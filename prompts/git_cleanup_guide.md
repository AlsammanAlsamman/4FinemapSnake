# Git Repository Cleanup Guide

## Problem
Repository has too many files tracked, causing slow Git operations and large repository size.

## Solution Steps

### 1. Update .gitignore
Create a restrictive `.gitignore` that ignores everything except essential files:

```gitignore
# Ignore everything by default
*

# Keep essential files
!.gitignore
!README.md
!steps.sh
!submit.sh
!test_loci_detection.sh

# Keep essential directories
!utils/
!utils/**
!scripts/
!scripts/**
!rules/
!rules/**
!configs/
!configs/**

# Ignore Python cache
__pycache__/
*.pyc
*.pyo

# Ignore large directories
results/
logs/
inputs/
outputs/
.snakemake/
```

### 2. Remove Files from Git Tracking
```bash
git rm -r --cached .
git add .gitignore steps.sh submit.sh test_loci_detection.sh utils/ scripts/ rules/ configs/ README.md
git commit -m "Keep only essential files"
```

### 3. Create Fresh Repository (if needed)
If history is too large, start fresh:

```bash
# Backup current .git folder (optional)
# Remove old git history
Remove-Item -Force -Recurse .git

# Initialize new repository
git init
git add .
git commit -m "Initial commit with essential files only"

# Add remote and force push
git remote add origin https://github.com/USERNAME/REPO.git
git branch -M main
git push -u origin main --force
```

## Result
- Only essential files tracked
- Fast Git operations
- Small repository size
- Clean history
