"""Full validation suite (47 checks, includes the parameter fits; ~1–2 min)."""
from urea_nitrate_cell.validation import run_all


def test_all_validation_checks_pass():
    rows = run_all()
    failed = [r["name"] for r in rows if not r["ok"]]
    assert len(rows) == 47
    assert not failed, failed
