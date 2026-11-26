"""Run the worked examples with hidden figures during validation."""

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def run_driver_changed(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.delenv("P2_OUTPUT_DIR", raising=False)

    def run(folder, filename, changes, names=()):
        path = ROOT / folder / filename
        text = path.read_text()
        for old, new in changes:
            assert text.count(old) == 1, old
            text = text.replace(old, new)
        namespace = {"__file__": str(path), "__name__": "__main__"}
        try:
            exec(compile(text, str(path), "exec"), namespace)
        except ValueError as error:
            return None, str(error)
        finally:
            plt.close("all")
        return {name: getattr(namespace[name], "tolist", lambda: namespace[name])()
                for name in names}, None

    return run


@pytest.fixture
def run_driver(run_driver_changed):
    def run(folder, filename, names):
        values, error = run_driver_changed(folder, filename, (), names)
        assert error is None, error
        return values, None

    return run
