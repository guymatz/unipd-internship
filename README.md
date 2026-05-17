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
