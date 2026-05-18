# Setup
## Setup Python

### "Production"
`pip install -r requirements.txt`

### Development
`pip install -r requirements-dev.txt`

#### For Jupyter Development
`pip install -r requirements-jupyter.txt`

## Setup pre-commit
```
brew install pre-commit gitleaks
pre-commit install
```
### Validate pre-commit installation
`pre-commit run --verbose --all-files`

## Other git Setup

### Add tests submodule, if not already done and if necessary
git submodule add git@github.com:guymatz/unipd-internship-tests.git tests
git config set submodule.recurse true
git config set submodule.propagateBranches true

### Add a post-checkout hook to keep tests submodule in sync
ln -s -f git-hooks/post-checkout .git/hooks/post-checkout
