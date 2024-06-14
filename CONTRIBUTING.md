# How to contribute to MulensModel?

Any contribution to the MulensModel package is welcomed, even if it is just another test or example on how to use MulensModel.
We note contributions in [AUTHORS.md](AUTHORS.md) file.

### Have you found a bug?

If you see a bug in the code or would like to see a new feature, then please [open a new issue](https://github.com/rpoleski/MulensModel/issues/new). If you don't think an issue is the right category, then start [a new discussion](https://github.com/rpoleski/MulensModel/discussions/new).

### Would you like to see a new feature?

We're also open to suggestions, new ideas, and requests. If you have them, then just [open a new issue](https://github.com/rpoleski/MulensModel/issues/new). We have a long list of changes that we would like to make and your feedback may help us prioritize them. A fresh look at the code is always useful!

### Do you have code to contribute?

Having a code that you want to share is great! We follow coding conventions (e.g., [PEP8](https://github.com/rpoleski/MulensModel/issues/new), [sphinx](https://www.sphinx-doc.org/en/master/)-compatible docstrings for public functions, [camel case](https://en.wikipedia.org/wiki/Snake_case) names of classes, [snake case](https://en.wikipedia.org/wiki/Snake_case) names of attributes and methods, etc.) but you don't have to follow all of them right away. Just open a new [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) and we will have some discussion on changes that need to be done before merging. Please note that we expect that you will provide a unit test for every new public function that is not plotting. If the function is plotting, then we expect a new example. Sometimes there is significant discussion befor the pull request is accepted (see some [closed ones](https://github.com/rpoleski/MulensModel/pulls?q=is%3Apr+is%3Aclosed) for examples).

An important aspect of adding new code is writing tests that can be run in automated way. The tests are in `source/MulensModel/tests/` directory and you can run them by invoking `pytest` command (it requires pytest package to be installed). If you don't know how to write tests for your code, don't worry - we will help you. We expect that new code is tested in nearly 100% when it's merged.

### Weekly hack sessions

As of Nov 2022, we have regular hack sessions (just 2h each week) to work on MulensModel. Contact us if you are interested in joining.


Thanks!

Radek & Jennifer

