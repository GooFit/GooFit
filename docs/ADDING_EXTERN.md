# Adding external dependencies

The CMake build system and git make adding dependencies easy. The following methods can be used.

## New submodule

Most of the smaller libraries are added using submodules. To add a library from GitHub, with user `username` and package name `package`:
```bash
git submodule add ../../username/package.git extern/package
```

Go into the newly downloaded package in `/extern` and checkout the approriate tag/commit. Remember: when you are in a submodule, git commands affect that submodule as if it were a normal git repository (even up to commit and push), while if you are above that submodule in the main git folder, git commands affect the main repository. So you will want to do `git add extern/package` to add the correct commit.

You will also need to modify the `CMakeLists.txt` file. For a header only library example:
```cmake
include_directories("${PROJECT_SOURCE_DIR}/extern/project/include")
```

For a CMake library, use `add_subdirectory(extern/package)` instead, and then link to the targets provided.

## DownloadPackage

The DownloadPackage cmake file is part of the `/cmake` submodule. The `AddX.cmake` files are examples of how to use this method.

## FindPackage

This is the "standard" method for adding cmake packages, but requires that the library already be available. It is used for ROOT and CUDA currently.


