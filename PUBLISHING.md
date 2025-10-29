# Publishing PyBFS to PyPI

This guide explains how to publish the PyBFS package to PyPI.

## Prerequisites

1. **PyPI Account**: Create an account at https://pypi.org/account/register/
2. **Test PyPI Account** (recommended for testing): Create an account at https://test.pypi.org/account/register/
3. **Install build tools**:
   ```bash
   pip install build twine
   ```

## Steps to Publish

### 1. Update Version Number

Before publishing, update the version in:
- `pybfs/__init__.py` (update `__version__`)
- `pyproject.toml` (update `version = ...`)

Use [Semantic Versioning](https://semver.org/):
- **Major version** (1.0.0): Breaking changes
- **Minor version** (0.2.0): New features, backward compatible
- **Patch version** (0.1.1): Bug fixes, backward compatible

### 2. Build Distribution Packages

```bash
# Clean any previous builds
rm -rf dist/ build/ *.egg-info

# Build source and wheel distributions
python -m build
```

This creates:
- `dist/pybfs-0.1.0.tar.gz` (source distribution)
- `dist/pybfs-0.1.0-py3-none-any.whl` (wheel distribution)

### 3. Check Package (Optional but Recommended)

```bash
# Check the built package
twine check dist/*
```

### 4. Test on Test PyPI (Recommended)

Before publishing to production PyPI, test on Test PyPI:

```bash
# Upload to Test PyPI
twine upload --repository testpypi dist/*

# Test installation
pip install --index-url https://test.pypi.org/simple/ pybfs
```

You'll be prompted for your Test PyPI username and password.

### 5. Publish to Production PyPI

Once tested, publish to production PyPI:

```bash
twine upload dist/*
```

You'll be prompted for your PyPI username and password.

**Note**: For better security, use an API token instead of password:
1. Go to https://pypi.org/manage/account/
2. Create an API token under "API tokens"
3. Use the token as the username and `__token__` as the password (or use `twine upload --username __token__ --password <token>`)

### 6. Verify Installation

After publishing, verify the package can be installed:

```bash
pip install pybfs
```

### 7. Update GitHub (Recommended)

Tag the release in git:

```bash
git add .
git commit -m "Release version 0.1.0"
git tag v0.1.0
git push origin main --tags
```

## Updating the Package

To release a new version:

1. Update version numbers (see Step 1)
2. Update `CHANGELOG.md` if you maintain one
3. Follow steps 2-7 above

## Common Issues

### Package name already taken
- If "pybfs" is taken, consider "py-bfs" or "pybfs-hydro" (and update `pyproject.toml`)

### Authentication errors
- Ensure you're using the correct PyPI credentials
- Consider using API tokens instead of passwords

### Upload errors
- Check that you've incremented the version number
- Verify `pyproject.toml` is correctly formatted

## Resources

- [PyPI User Guide](https://packaging.python.org/en/latest/guides/distributing-packages-using-setuptools/)
- [Python Packaging Guide](https://packaging.python.org/)
- [Semantic Versioning](https://semver.org/)

