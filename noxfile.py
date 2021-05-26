import nox
import multiprocessing

@nox.session
def test(session):
    session.env["CMAKE_BUILD_PARALLEL_LEVEL"] = str(multiprocessing.cpu_count())
    session.install("-v", ".[dev]", silent=False)
    session.run("pytest")
