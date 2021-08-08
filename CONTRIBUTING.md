# Contributing to Celest

## Getting Started

Before you begin:
* Have you read the [code of conduct](./CODE_OF_CONDUCT.md)?
* See if an [issue exists](https://github.com/JaiWillems/Celest/issues) for the change you want to make.
* See if we are [accepting controbutions](#Types-of-Contributions) for your type of issue.

### Don't see your issue? Open one

If you have a change that you would like to make, open an issue using a [template](https://github.com/JaiWillems/Celest/issues/new/choose). This will be the place to discuss the problem you want to fix.

### Ready to make a change? Fork the repo

Fork using GitHub Desktop:
* [Getting started with GitHub Desktop](https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/overview/getting-started-with-github-desktop) will guide you through setting up Desktop.
* Once Desktop is set up, you can use it to [fork the repo](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop)!

Fork using the command line:
* [Fork the repo](https://docs.github.com/en/get-started/quickstart/fork-a-repo#fork-an-example-repository) so that you can make changes without affecting the original project until you are ready to merge them.

Fork with [GitHub Codespaces](https://github.com/features/codespaces):
* [Fork, edit, and preview](https://docs.github.com/en/codespaces/developing-in-codespaces/creating-a-codespace) using [GitHub Codespaces](https://github.com/features/codespaces) without having to install and run the project locally.

### Make your update:

Make the changes to the file(s) you'd like to update.

### Open a pull request

When you're done making changes and you'd like to propose them for review, use the [pull request template](#Pull-request-template) to open your pull request.

### Submit your pull request & get it reviewed
* Once you submit your pull request, others from the Celest community will review it with you. The first thing you're going to want to do is a self-review.
* After that, we may have questions, check back on your pull request to keep up with the conversation.
* Did you have an issue, like a merge conflict? Check out the [git tutorial](https://lab.github.com/githubtraining/managing-merge-conflicts) on how to resolve merge conflicts and other issues.

### Your pull request is merged!
Congratulations! The Celest community thanks you!

Once your pull request is merged, you will proudly be listed as a contributor in the [contributor chart](https://github.com/JaiWillems/Celest/graphs/contributors).

## Types of Contributions

You can contribute to Celest in several ways.

### Issues

Issues are used to track tasks that contributors can help with. If an issue hasn't received a reply from a Celest maintainer, then we haven't reviewed the issue yet and you shouldn't begin to work on it.

If you have found something in the code or documentation that should be updated, search open issues to see if someone else has reported the same thing. If it's something new, open an issue using a [template](https://github.com/JaiWillems/Celest/issues/new/choose). We'll use the issue to have a conversation about the problem you want to fix.

### Pull Requests

A pull request is a way to suggest changes in our repository.

When we merge those changes to the current develop branch, they will be deployed in the upcoming Celest release. To learn more about opening a pull request in this repo, see [Opening a Pull Request](#Opening-a-Pull-Request) below.

### Support

We are a small team working hard to keep up with the demands of a continuously changing product. Unfortunately, we just can't help with support questions in this repository. If you are experiencing a problem with GitHub, unrelated to Celest, please [contact GitHub Support directly](https://support.github.com/request). Any issues, discussions, or pull requests opened here requesting support will be given information about how to contact GitHub Support, then closed and locked.

If you're having trouble with your GitHub account, contact [Support](https://support.github.com/contact?tags=docs-contributing-guide).

## Starting With an Issue

You can browse existing issues to find something that needs help!

### Labels
Labels can help you find an issue you'd like to help with.

* The [bug label](https://github.com/JaiWillems/Celest/labels/bug) is for problems where something isn't working.
* The [documentation label](https://github.com/JaiWillems/Celest/labels/documentation) is for improvements or additions to documentation.
* The [duplicate label](https://github.com/JaiWillems/Celest/labels/duplicate) is for any issue or pull request that already exists.
* The [enhancement label](https://github.com/JaiWillems/Celest/labels/enhancement) is for a new feature or request.
* The [good first issue label](https://github.com/JaiWillems/Celest/labels/good%20first%20issue) is a good issue for first-time contributors to work on.
* The [help wanted label](https://github.com/JaiWillems/Celest/labels/help%20wanted) means extra attention is needed to the issue.
* The [invalid label](https://github.com/JaiWillems/Celest/labels/invalid) is for issues where there appear to be invalid results.
* The [question label](https://github.com/JaiWillems/Celest/labels/question) is to request further information.
* The [wontfix](https://github.com/JaiWillems/Celest/labels/wontfix) label indicates the issue will not be worked on.

## Opening a Pull Request

You can use the GitHub user interface for small changes like fixing a typo or updating a readme. You can also fork the repo and then clone it locally, to view changes and run your tests on your machine.

## Reviewing

We review every single pull request. The purpose of reviews is to create the best content we can for people who use Celest.

* Reviews are always respectful, acknowledging that everyone did the best possible job with the knowledge that they had at the time.
* Reviews discuss content, not the person who created it.
* Reviews are constructive and start conversations around feedback.

### Self Review

You should always review your pull request first.

For code changes, make sure that you:
* Confirm that the changes meet the user experience and goals outlined in the content design plan (if there is one).
* Ensure that technical changes are validated against a truth. Proof of validation for non-trivial outputs will need to be presented before a pull request is approved for merging.
* Copy-edit the added documentation for grammar, spelling, and adherence to the style guide.
* If there are any failing checks in your pull request, troubleshoot them until they're all passing.

For all code additions, it is encouraged to:
* [Type annotate](https://docs.python.org/3/library/typing.html) new function and methods.
* Ensure code meets the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/). A useful online checker that will point out a number of PEP8 infractions can be found [here](http://pep8online.com/).
* Ensure documentation strings are created or updated per [Numpy's style guide](https://numpydoc.readthedocs.io/en/latest/format.html).

### Pull Request Template
When ready to open a pull request on your issue, fill out the pull request template to ensure adequate information is incorporated for faster code reviews.

### Suggested Changes

We may ask for changes to be made before a pull request can be merged, either using pull request comments. You can apply suggested changes directly through the UI. You can make any other changes in your fork, then commit them to your branch.

As you update your pull request and apply changes, mark each conversation as resolved.