"""Mascope data-curation engine — the WRITE side of the API.

`io_mascope` is the read boundary (peaks, matches, listings). This module is its
mutating counterpart: create / rename / copy / move / delete of **workspaces,
datasets, batches and sample items**, plus reagent-match re-runs. None of these
endpoints are in `mascope-sdk` (which is read-only + `match_compounds`); they were
reverse-engineered from the Mascope web app's own API calls and validated live
(see `docs/DATA_CURATION.md` for the endpoint spec + how it was derived).

Design, mirroring `io_mascope`:

  * ONE HTTP boundary. Everything goes through ``CurationClient._request`` -> the
    same ``{url}/api/{path}`` shape, bearer auth, Cloudflare-WAF retry, and the
    ``{"data": ...}`` envelope unwrap the read side uses.
  * SAFE BY DEFAULT. Every mutating call is recorded in ``client.actions``. A
    client built with ``dry_run=True`` PLANS every mutation (records + returns the
    request it *would* send) and sends nothing — so you can preview a whole
    reorganisation, show it to the user, then re-run for real. Destructive deletes
    additionally require an explicit ``confirm=True`` (or dry-run) — a bare
    ``delete_*`` refuses to fire.
  * NAME-DRIVEN. ``resolve_*_id`` map human names -> ids (exact-id > exact-name >
    substring, raising on 0/ambiguous, like the SDK), and ``ensure_*`` are
    idempotent get-or-create helpers — so Claude can drive the whole thing by the
    names a scientist actually uses ("move the zero-air samples into a Blanks
    batch") without hand-managing opaque ids.

The engine is deliberately thin and composable: it exposes the primitives, and the
sorting/curation *policy* (which batch a sample belongs in, how to split a campaign)
lives in the calling session, not here.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import pandas as pd
import requests

from peaky.io import io_mascope as IO

__version__ = "0.1.0"

# Methods that change server state. dry-run plans these; reads always execute.
_MUTATING = frozenset({"POST", "PATCH", "PUT", "DELETE"})
# Mascope polarity tokens (batch create field is `sample_batch_polarity`).
POLARITIES = ("+", "-", "+-")


class CurationError(RuntimeError):
    """A curation request failed, or a guard (confirm/polarity) rejected a call."""


@dataclass
class Action:
    """One recorded mutation — the audit trail + the dry-run plan.

    ``executed`` is False for a planned (dry-run or unconfirmed) action; the
    ``result`` holds the server's response data once it actually ran.
    """
    op: str                       # high-level verb, e.g. 'create_batch'
    method: str                   # HTTP method
    path: str                     # api path (no /api/ prefix)
    body: dict | None = None      # request JSON body
    params: dict | None = None    # query string
    executed: bool = False
    result: object = None
    note: str = ""

    def describe(self) -> str:
        tag = "DONE" if self.executed else "PLAN"
        extra = f"  {self.note}" if self.note else ""
        return f"[{tag}] {self.op}: {self.method} /{self.path}{extra}"


@dataclass
class CurationClient:
    """Authenticated write client for one Mascope server.

    Build it with :meth:`from_env` (reads ``MASCOPE_URL`` / ``MASCOPE_ACCESS_TOKEN``
    from the same ``.env`` search order as the read side). Pass ``dry_run=True`` to
    PLAN every mutation without sending it.
    """
    url: str
    token: str
    verify_ssl: bool = False
    dry_run: bool = False
    service_name: str = "mascope_sdk"
    timeout: tuple = (15, 120)
    cookie: str | None = None       # session cookie for WRITE endpoints (see below)
    actions: list[Action] = field(default_factory=list)

    # NOTE on AUTH (validated 2026-07-15): the bearer API token authorizes the
    # READ endpoints + match_compounds, but the Mascope server gates the
    # WRITE/curation endpoints (create/copy/move/rename/delete) behind an
    # interactive SESSION COOKIE — the web app sends requests with
    # ``withCredentials:true`` and no Bearer at all. A token-only write returns
    # 401 "Please sign in to the Mascope". So to actually mutate, pass ``cookie=``
    # (or set ``MASCOPE_SESSION_COOKIE``): the full ``Cookie:`` header value copied
    # from a logged-in browser session (DevTools → Network → any /api request →
    # Request Headers → Cookie). The bearer is still sent alongside (harmless).

    # -- construction --------------------------------------------------------
    @classmethod
    def from_env(cls, env_path: str | None = None, *, dry_run: bool = False,
                 verify_ssl: bool = False, cookie: str | None = None
                 ) -> "CurationClient":
        """Resolve creds via the shared ``io_mascope`` ``.env`` search order (or the
        process env). The optional session ``cookie`` (needed for writes) comes from
        the argument or ``$MASCOPE_SESSION_COOKIE``. Raises with an actionable
        message when no URL/token is found."""
        from dotenv import load_dotenv
        path = IO._find_env(env_path)
        load_dotenv(path)
        url = os.environ.get("MASCOPE_URL")
        tok = os.environ.get("MASCOPE_ACCESS_TOKEN")
        if not url or not tok:
            raise CurationError(
                f"MASCOPE_URL / MASCOPE_ACCESS_TOKEN not found (looked in {path}). "
                "Fill the repo-root .env (or ~/.mascope/.env), set $MASCOPE_ENV, or "
                "export the two variables.")
        return cls(url=url.rstrip("/"), token=tok, verify_ssl=verify_ssl,
                   dry_run=dry_run,
                   cookie=cookie or os.environ.get("MASCOPE_SESSION_COOKIE") or None)

    # -- the one HTTP boundary ----------------------------------------------
    def _request(self, method: str, path: str, *, json: dict | None = None,
                 params: dict | None = None, op: str = "", note: str = ""):
        """Send one request (or PLAN it, if dry-run + mutating). Returns the
        unwrapped ``data`` payload (``None`` on empty/204). Mutating calls are
        appended to ``self.actions`` either way."""
        method = method.upper()
        mutating = method in _MUTATING
        act = Action(op=op or method.lower(), method=method, path=path,
                     body=json, params=params, note=note)

        if mutating and self.dry_run:
            act.executed = False
            act.note = (note + " (dry-run: not sent)").strip()
            self.actions.append(act)
            return {"_planned": True, "method": method, "path": path,
                    "body": json, "params": params}

        full = f"{self.url}/api/{path}"
        headers = {"Authorization": f"Bearer {self.token}",
                   "X-Service-Name": self.service_name}
        if self.cookie:                 # session auth — required by write endpoints
            headers["Cookie"] = self.cookie
        if mutating and not self.cookie:
            raise CurationError(
                "This is a WRITE endpoint; the Mascope server requires a logged-in "
                "session cookie for it (the bearer token only authorizes reads + "
                "match_compounds). Pass cookie=... / set $MASCOPE_SESSION_COOKIE with "
                "the Cookie header from a logged-in browser session, or build the "
                "client with dry_run=True to preview without sending.")

        def _send():
            r = requests.request(method, full, headers=headers, json=json,
                                 params=params, timeout=self.timeout,
                                 verify=self.verify_ssl)
            self._raise_for_status(r)
            return r

        r = IO._with_waf_retry(_send)
        data = None
        if r.content:
            try:
                data = r.json().get("data")
            except ValueError:
                data = None
        if mutating:
            act.executed = True
            act.result = data
            self.actions.append(act)
        return data

    @staticmethod
    def _raise_for_status(r: requests.Response) -> None:
        """Raise ``CurationError`` with the server's own message on a non-2xx."""
        if r.ok:
            return
        msg = ""
        try:
            body = r.json()
            if isinstance(body, dict):
                det = body.get("detail")
                if isinstance(det, dict):
                    msg = det.get("error_message") or str(det)
                elif det:
                    msg = str(det)
                else:
                    msg = body.get("message") or body.get("error") or ""
        except ValueError:
            msg = (r.text or "")[:200]
        raise CurationError(f"HTTP {r.status_code} on {r.request.method} "
                            f"{r.url}: {msg or 'no message'}")

    # -- convenience verbs ---------------------------------------------------
    def _get(self, path, **params):
        return self._request("GET", path, params=params or None)

    # =======================================================================
    # READ helpers (name -> id resolution needs them; also handy standalone)
    # =======================================================================
    def list_workspaces(self) -> pd.DataFrame:
        return pd.DataFrame(self._get("workspaces") or [])

    def list_datasets(self, workspace_id: str) -> pd.DataFrame:
        return pd.DataFrame(self._get(f"workspaces/{workspace_id}/datasets") or [])

    def list_batches(self, dataset_id: str) -> pd.DataFrame:
        return pd.DataFrame(self._get("sample/batches", dataset_id=dataset_id) or [])

    def list_samples(self, sample_batch_id: str) -> pd.DataFrame:
        return pd.DataFrame(self._get("samples", sample_batch_id=sample_batch_id) or [])

    # -- name -> id ----------------------------------------------------------
    @staticmethod
    def _resolve(df: pd.DataFrame, name_or_id: str, *, id_col: str, name_col: str,
                 kind: str) -> str:
        """exact-id > exact-name (case-insensitive) > substring; raise on 0/>1
        (mirrors the SDK's resolution + loud ambiguity)."""
        if df is None or not len(df):
            raise CurationError(f"no {kind}s on the server to match {name_or_id!r}")
        key = str(name_or_id)
        if id_col in df.columns:
            hit = df[df[id_col].astype(str) == key]
            if len(hit) == 1:
                return str(hit.iloc[0][id_col])
        names = df[name_col].astype(str)
        hit = df[names.str.casefold() == key.casefold()]
        if not len(hit):
            hit = df[names.str.casefold().str.contains(key.casefold(), regex=False)]
        if len(hit) == 1:
            return str(hit.iloc[0][id_col])
        if not len(hit):
            avail = ", ".join(sorted(names)[:30])
            raise CurationError(f"no {kind} matching {name_or_id!r}. Available: {avail}")
        opts = ", ".join(hit[name_col].astype(str).tolist())
        raise CurationError(f"{len(hit)} {kind}s match {name_or_id!r}; be specific: {opts}")

    def resolve_workspace_id(self, name_or_id: str) -> str:
        return self._resolve(self.list_workspaces(), name_or_id,
                             id_col="workspace_id", name_col="workspace_name",
                             kind="workspace")

    def resolve_dataset_id(self, workspace: str, name_or_id: str) -> str:
        wid = self.resolve_workspace_id(workspace)
        return self._resolve(self.list_datasets(wid), name_or_id,
                             id_col="dataset_id", name_col="dataset_name",
                             kind="dataset")

    def resolve_batch_id(self, dataset_id: str, name_or_id: str) -> str:
        return self._resolve(self.list_batches(dataset_id), name_or_id,
                             id_col="sample_batch_id", name_col="sample_batch_name",
                             kind="batch")

    # =======================================================================
    # WORKSPACES
    # =======================================================================
    def create_workspace(self, name: str, description: str = "") -> dict:
        return self._request("POST", "workspaces", op="create_workspace",
                             json={"workspace_name": name,
                                   "workspace_description": description},
                             note=f"name={name!r}")

    def update_workspace(self, workspace_id: str, *, name: str | None = None,
                         description: str | None = None) -> dict:
        body: dict = {"workspace_id": workspace_id}
        if name is not None:
            body["workspace_name"] = name
        if description is not None:
            body["workspace_description"] = description
        return self._request("PATCH", f"workspaces/{workspace_id}",
                             op="update_workspace", json=body, note=f"name={name!r}")

    def delete_workspace(self, workspace_id: str, *, confirm: bool = False) -> dict:
        self._guard_delete("workspace", workspace_id, confirm)
        return self._request("DELETE", f"workspaces/{workspace_id}",
                             op="delete_workspace", note=f"id={workspace_id}")

    # =======================================================================
    # DATASETS
    # =======================================================================
    def create_dataset(self, workspace_id: str, name: str,
                       description: str = "") -> dict:
        return self._request("POST", f"workspaces/{workspace_id}/datasets",
                             op="create_dataset",
                             json={"dataset_name": name,
                                   "dataset_description": description},
                             note=f"name={name!r} ws={workspace_id}")

    def update_dataset(self, workspace_id: str, dataset_id: str, *,
                       name: str | None = None, description: str | None = None,
                       dataset_type: str | None = None) -> dict:
        body: dict = {"dataset_id": dataset_id}
        if name is not None:
            body["dataset_name"] = name
        if description is not None:
            body["dataset_description"] = description
        if dataset_type is not None:
            body["dataset_type"] = dataset_type
        return self._request("PATCH", f"workspaces/{workspace_id}/datasets/{dataset_id}",
                             op="update_dataset", json=body, note=f"name={name!r}")

    def delete_dataset(self, workspace_id: str, dataset_id: str, *,
                       confirm: bool = False) -> dict:
        self._guard_delete("dataset", dataset_id, confirm)
        return self._request("DELETE",
                             f"workspaces/{workspace_id}/datasets/{dataset_id}",
                             op="delete_dataset", note=f"id={dataset_id}")

    def move_dataset(self, dataset_id: str, source_workspace_id: str,
                     target_workspace_id: str) -> dict:
        """Move a dataset from one workspace to another (POST .../datasets/{id}/move)."""
        return self._request(
            "POST", f"workspaces/{source_workspace_id}/datasets/{dataset_id}/move",
            op="move_dataset",
            json={"source_workspace_id": source_workspace_id,
                  "target_workspace_id": target_workspace_id},
            note=f"{source_workspace_id} -> {target_workspace_id}")

    # =======================================================================
    # BATCHES
    # =======================================================================
    def create_batch(self, dataset_id: str, name: str, polarity: str,
                     description: str = "",
                     target_collection_ids: list | None = None) -> dict:
        """Create an empty sample batch in a dataset. `polarity` in {'+','-','+-'}."""
        if polarity not in POLARITIES:
            raise CurationError(f"polarity must be one of {POLARITIES}, got {polarity!r}")
        return self._request(
            "POST", "sample/batches", op="create_batch",
            json={"sample_batch_name": name, "sample_batch_description": description,
                  "sample_batch_polarity": polarity, "dataset_id": dataset_id,
                  "target_collection_ids": list(target_collection_ids or [])},
            note=f"name={name!r} pol={polarity} ds={dataset_id}")

    def update_batch(self, sample_batch_id: str, *, name: str | None = None,
                     description: str | None = None) -> dict:
        body: dict = {"sample_batch_id": sample_batch_id}
        if name is not None:
            body["sample_batch_name"] = name
        if description is not None:
            body["sample_batch_description"] = description
        return self._request("PATCH", f"sample/batches/{sample_batch_id}",
                             op="update_batch", json=body, note=f"name={name!r}")

    def delete_batch(self, sample_batch_id: str, *, confirm: bool = False) -> dict:
        self._guard_delete("batch", sample_batch_id, confirm)
        return self._request("DELETE", f"sample/batches/{sample_batch_id}",
                             op="delete_batch", note=f"id={sample_batch_id}")

    def copy_batch(self, sample_batch_id: str, target_dataset_id: str, *,
                   name: str, description: str = "") -> dict:
        """Copy a whole batch (with its samples) into a dataset under a new name."""
        return self._request(
            "POST", f"sample/batches/{sample_batch_id}/copy", op="copy_batch",
            json={"dataset_id": target_dataset_id, "sample_batch_name": name,
                  "sample_batch_description": description},
            note=f"-> ds={target_dataset_id} as {name!r}")

    def import_samples(self, sample_batch_id: str, sample_items: list) -> dict:
        """Attach existing sample items to a batch (POST .../import)."""
        return self._request("POST", f"sample/batches/{sample_batch_id}/import",
                             op="import_samples", json={"sample_items": sample_items},
                             note=f"n={len(sample_items)}")

    def rematch_batch(self, sample_batch_id: str, *, full_remove: bool = False,
                      force: bool = False) -> dict:
        """Re-run the target-list match for a batch (POST /match/rematch/batch/{id})."""
        return self._request("POST", f"match/rematch/batch/{sample_batch_id}",
                             op="rematch_batch", json={},
                             params={"full_remove": full_remove, "force": force},
                             note=f"full_remove={full_remove} force={force}")

    # =======================================================================
    # SAMPLE ITEMS (the measurement rows) — the workhorse for sorting data
    # =======================================================================
    def update_sample(self, sample_item_id: str, **fields) -> dict:
        """Patch a sample's metadata (name, type, attributes, ...). The server
        wraps the payload as ``{"sample_item": {...}}``."""
        if not fields:
            raise CurationError("update_sample needs at least one field to change")
        return self._request("PATCH", f"sample/items/{sample_item_id}",
                             op="update_sample", json={"sample_item": dict(fields)},
                             note=f"fields={sorted(fields)}")

    def copy_samples(self, sample_item_ids: list, target_batch_id: str) -> dict:
        """Copy sample items INTO another batch (originals stay put)."""
        return self._request("POST", "sample/items/copy", op="copy_samples",
                             json={"sample_batch_id": target_batch_id,
                                   "sample_item_ids": list(sample_item_ids)},
                             note=f"n={len(sample_item_ids)} -> {target_batch_id}")

    def move_samples(self, sample_item_ids: list, target_batch_id: str) -> dict:
        """Move sample items into another batch (removed from the source)."""
        return self._request("POST", "sample/items/move", op="move_samples",
                             json={"sample_batch_id": target_batch_id,
                                   "sample_item_ids": list(sample_item_ids)},
                             note=f"n={len(sample_item_ids)} -> {target_batch_id}")

    def delete_samples(self, sample_item_ids: list, *, confirm: bool = False) -> dict:
        self._guard_delete("sample", f"{len(sample_item_ids)} items", confirm)
        return self._request("POST", "sample/items/delete", op="delete_samples",
                             json={"sample_item_ids": list(sample_item_ids)},
                             note=f"n={len(sample_item_ids)}")

    # =======================================================================
    # RAW FILES (ingestion side)
    # =======================================================================
    def delete_files(self, sample_file_ids: list, *, confirm: bool = False) -> dict:
        self._guard_delete("file", f"{len(sample_file_ids)} files", confirm)
        return self._request("POST", "sample/files/delete", op="delete_files",
                             json={"sample_file_ids": list(sample_file_ids)},
                             note=f"n={len(sample_file_ids)}")

    def reprocess_files(self, sample_file_ids: list) -> dict:
        return self._request("POST", "sample/files/reprocess", op="reprocess_files",
                             json={"sample_file_ids": list(sample_file_ids)},
                             note=f"n={len(sample_file_ids)}")

    # =======================================================================
    # Idempotent get-or-create (name-driven curation)
    # =======================================================================
    def ensure_workspace(self, name: str, description: str = "") -> str:
        """Return the id of the workspace named `name`, creating it if absent.
        In dry-run when creation is needed, returns a planning sentinel string."""
        try:
            return self.resolve_workspace_id(name)
        except CurationError:
            res = self.create_workspace(name, description)
            return _new_id(res, "workspace_id")

    def ensure_dataset(self, workspace: str, name: str, description: str = "") -> str:
        wid = self.resolve_workspace_id(workspace)
        try:
            return self._resolve(self.list_datasets(wid), name, id_col="dataset_id",
                                name_col="dataset_name", kind="dataset")
        except CurationError:
            res = self.create_dataset(wid, name, description)
            return _new_id(res, "dataset_id")

    def ensure_batch(self, dataset_id: str, name: str, polarity: str,
                     description: str = "") -> str:
        try:
            return self.resolve_batch_id(dataset_id, name)
        except CurationError:
            res = self.create_batch(dataset_id, name, polarity, description)
            return _new_id(res, "sample_batch_id")

    # -- plan / audit --------------------------------------------------------
    def plan(self) -> list[Action]:
        """The recorded mutations (dry-run = what WOULD happen; live = what did)."""
        return list(self.actions)

    def summary(self) -> str:
        if not self.actions:
            return "(no mutations recorded)"
        return "\n".join(a.describe() for a in self.actions)

    def _guard_delete(self, kind: str, ident: str, confirm: bool) -> None:
        """A destructive delete must be explicitly confirmed OR planned (dry-run).
        Blocks an accidental unconfirmed delete even in live mode."""
        if self.dry_run or confirm:
            return
        raise CurationError(
            f"refusing to delete {kind} {ident} without confirm=True "
            "(or build the client with dry_run=True to preview it).")


def _new_id(result, id_col: str) -> str:
    """Pull the new id out of a create response, tolerating the dry-run sentinel."""
    if isinstance(result, dict):
        if result.get("_planned"):
            return f"<planned:{id_col}>"
        if id_col in result:
            return str(result[id_col])
    raise CurationError(f"create response lacked {id_col}: {result!r}")
