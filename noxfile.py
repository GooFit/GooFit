import multiprocessing
import shutil
from pathlib import Path

import nox


@nox.session
def test(session):
    session.env["CMAKE_BUILD_PARALLEL_LEVEL"] = str(multiprocessing.cpu_count())
    if Path("_skbuild").exists():
        shutil.rmtree("_skbuild")
    session.install("-v", ".[dev]", silent=False)
    session.run("pytest")
