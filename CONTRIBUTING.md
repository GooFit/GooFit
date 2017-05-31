# How to contribute to GooFit?

## Adding a library

See [adding a package](doc/ADDING_EXTERN.md).

## C++

You should verify that the LLVM 4.0 based requirements pass. To fix the clang-tidy requirements:

```bash
cmake -DGOOFIT_TIDY_FIX=ON ..
make -j1
```

To run `clang-format` to correct source code spacing:

```bash
git ls-files -- '*.cu' '*.cc' '*.h' '*.cpp' | xargs clang-format -i -style=file
```

The modernize script should also not make any changes:

```bash
git ls-files -- '*.cu' '*.cc' '*.h' '*.cpp' | xargs ModernizeGoofit.py
```

## Git

GooFit is hosted as a [git](http://git-scm.com) repository on [GitHub](https://github.com).

If you know git and GitHub, contributing to GooFit is easy:
fork and make pull requests.
Use the [git help](http://git-scm.com/documentation) and [GitHub help](https://help.github.com) pages
or Google if you don't remember how to do something or get an error.

If you want to make a small change you can use the GitHub to edit or create files on the GitHub web page as described [here](https://help.github.com/articles/creating-and-editing-files-in-your-repository). The advantage is that you don't need to know git and it's quick, the disadvantage is that you can't test your changes on your local computer, so this is not useful to code on GooFit. The built in Travis CI tests do not check CUDA code, so it is recomended that you test your code with
CUDA locally as well if you are making changes to computations in GooFit.

There are many great git and github tutorials out there, e.g.:
* http://gitimmersion.com
* http://try.github.io/

Git is very flexible and allows different workflows, we use the one used by most scientific open-source projects as described in detail e.g. by the [astropy project](http://docs.astropy.org/en/latest/development/workflow/development_workflow.html).

There's no point in repeating all that info specifically for GooFit here, we do hope the references above are enough to make it easy for you to learn how to contribute.

If you don't have time to learn about git and pull requests, we still do want your contribution to GooFit. In that case please make an issue on GitHub and link to your updated files or patches e.g. in [gists](https://gist.github.com).
Note that this way you won't get credit for you contribution though in the commit history, so we prefer you make a GitHub pull request or ask for help in a GitHub issue on where you got stuck (fork -> clone -> branch -> edit -> commit -> push -> pull request) with making the pull request.


## Contact

Please use [Issues](https://github.com/GooFit/GooFit/issues) for suggestions,
ask us questions on [Gitter](https://gitter.im/GooFit/Lobby),
and make [Pull Requests](https://github.com/GooFit/GooFit/pulls) when you are ready to contribute!
Please also let us know if you have an interesting analysis in GooFit so we can link to it
or mention it on our [public page](https://GooFit.github.io).
