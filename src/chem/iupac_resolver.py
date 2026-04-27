"""Multi-source IUPAC name resolution from validated SMILES."""

import json
import os
import socket
from urllib import parse, request


DEFAULT_TIMEOUT = 5
DEFAULT_STOUT_URL = "http://127.0.0.1:8000/predict"


def _empty_result(source):
    return {"source": source, "name": None, "success": False}


def _success_result(source, name):
    name = name.strip() if isinstance(name, str) else None
    if not name:
        return _empty_result(source)
    return {"source": source, "name": name, "success": True}


def _http_get_text(url, timeout=DEFAULT_TIMEOUT):
    req = request.Request(url, headers={"User-Agent": "ChemCanvas/1.0"})
    with request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace").strip()


def _http_post_json(url, payload, timeout=DEFAULT_TIMEOUT):
    data = json.dumps(payload).encode("utf-8")
    req = request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "ChemCanvas/1.0",
        },
        method="POST",
    )
    with request.urlopen(req, timeout=timeout) as response:
        text = response.read().decode("utf-8", errors="replace").strip()
    if not text:
        return {}
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return {"name": text}


def _rdkit_says_invalid(smiles):
    try:
        from rdkit import Chem, RDLogger
    except Exception:
        return False
    RDLogger.DisableLog("rdApp.error")
    try:
        return Chem.MolFromSmiles(smiles) is None
    finally:
        RDLogger.EnableLog("rdApp.error")


def resolve_with_cactus(smiles):
    try:
        encoded = parse.quote(smiles, safe="")
        url = "https://cactus.nci.nih.gov/chemical/structure/%s/iupac_name" % encoded
        name = _http_get_text(url)
        return _success_result("cactus", name)
    except Exception:
        return _empty_result("cactus")


def resolve_with_pubchem(smiles):
    try:
        import pubchempy as pcp

        previous_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(DEFAULT_TIMEOUT)
        try:
            compounds = pcp.get_compounds(smiles, namespace="smiles")
        finally:
            socket.setdefaulttimeout(previous_timeout)
        if not compounds:
            return _empty_result("pubchem")
        return _success_result("pubchem", compounds[0].iupac_name)
    except Exception:
        return _empty_result("pubchem")


def resolve_with_opsin(smiles):
    # OPSIN's public API is chemical-name to structure, not SMILES to IUPAC name.
    return _empty_result("opsin")


def resolve_with_stout(smiles):
    try:
        url = os.environ.get("STOUT_PREDICT_URL", DEFAULT_STOUT_URL)
        data = _http_post_json(url, {"smiles": smiles})
        for key in ("name", "iupac_name", "iupac", "prediction", "predicted_name"):
            if key in data:
                return _success_result("stout", data[key])
        return _empty_result("stout")
    except Exception:
        return _empty_result("stout")


def _dedupe_names(names):
    seen = set()
    deduped = []
    for name in names:
        key = name.strip().casefold()
        if key in seen:
            continue
        seen.add(key)
        deduped.append(name)
    return deduped


def resolve_iupac_name(smiles, use_stout=False):
    if not smiles or _rdkit_says_invalid(smiles):
        return {
            "input_smiles": smiles,
            "results": [],
            "valid_names": [],
            "preferred": None,
        }

    results = [
        resolve_with_cactus(smiles),
        resolve_with_pubchem(smiles),
        resolve_with_opsin(smiles),
    ]
    if use_stout:
        results.append(resolve_with_stout(smiles))

    valid = [result for result in results if result["success"] and result["name"]]
    valid_names = _dedupe_names([result["name"] for result in valid])

    return {
        "input_smiles": smiles,
        "results": results,
        "valid_names": valid_names,
        "preferred": valid_names[0] if valid_names else None,
    }


def get_iupac(smiles):
    return resolve_iupac_name(smiles)


def test_iupac_resolver():
    benzene = resolve_iupac_name("c1ccccc1")
    assert benzene["preferred"] is not None
    assert benzene["valid_names"]

    aspirin = resolve_iupac_name("CC(=O)Oc1ccccc1C(=O)O")
    names_by_source = {
        result["source"]: result["name"]
        for result in aspirin["results"]
        if result["success"] and result["name"]
    }
    assert "cactus" in names_by_source
    assert "pubchem" in names_by_source

    invalid = resolve_iupac_name("not a smiles")
    assert invalid["results"] == []
    assert invalid["valid_names"] == []
    assert invalid["preferred"] is None
