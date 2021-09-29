## Help setting up

See [`docs/SYSTEM_INSTALL.md`](../docs/SYSTEM_INSTALL.md) for help building and installing.

## Development

All formatting is done with [pre-commit][]. Install with brew (`brew install
pre-commit`) (macOS) or pipx/pip (any OS).

[pre-commit]: https://pre-commit.com


## Making a release (for developers)

* Make version bump PR
  * Bump the version numbers in setup.py and CMakeLists.txt
  * Set the `GOOFIT_TAG` to `"RELEASE"` in CMakeLists.txt
* Make GitHub release (web interface)
  * Github page -> Releases -> Draft a new release
  * Copy the release notes (should be up to date!) to the release (include link definitions)
  * CI should automatically append source package and upload docs / PyPI packages for you
* Return `GOOFIT_TAG` to `"dev"`
