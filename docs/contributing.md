# Contributing to PyBFS

Thank you for your interest in contributing to PyBFS!

## Understanding Git and GitHub Concepts

### The Big Picture

Contributing to PyBFS involves a collaborative workflow that protects the main codebase while giving you freedom to experiment. The workflow differs slightly depending on whether you're a core team member or an external contributor:

**For Core Team Members (with repository access):** You'll work directly with the main repository. Clone the PyBFS repository to your computer, then create a **branch** - a separate workspace where you can make changes without affecting the main code. As you work, you make **commits** to save snapshots of your progress. When you're done, you **push** your branch to the repository on GitHub and open a **pull request**, which is a formal proposal asking other team members to review and merge your changes into the main branch. Meanwhile, as other contributors' changes get merged, you'll periodically **pull** those updates from the main repository to stay synchronized.

**For External Contributors (without repository access):** You'll start by creating a **fork** (your personal copy) of the PyBFS repository on GitHub. You then **clone** that fork to your computer and work on branches within your fork. When you're done, you **push** your branch to your fork on GitHub and open a **pull request** to the main repository. You'll periodically **pull** updates from **upstream** (the original PyBFS repository) to stay synchronized with the project.

Both workflows ensure that everyone can contribute safely and collaboratively without stepping on each other's toes.

### Key Terms Explained

Let's break down each concept in detail:

### What is a Fork?

A **fork** is your personal copy of someone else's repository on GitHub. Think of it as getting your own playground version of the project where you can experiment and make changes without affecting the original. The original repository is often called the "upstream" repository.

**Why fork?** If you don't have direct write access to a repository, a fork gives you a space where you have full control to make changes, which you can later propose to be included in the original project.

**Note for PyBFS Team Members:** If you have admin or write access to the PyBFS repository, you typically don't need to fork - you can clone the main repository directly and work on branches.

### What is a Branch?

A **branch** is like a parallel version of your code. Imagine you're writing a book - the `main` branch is your published draft, but you create a new branch called `add-chapter-5` to work on a new chapter without disrupting the main draft.

**Why branch?** Branches let you work on multiple features simultaneously without them interfering with each other. If you're fixing a bug on one branch and adding a feature on another, they remain completely separate until you're ready to merge them.

### What is a Commit?

A **commit** is a snapshot of your changes at a specific point in time. Each commit has a message describing what changed. Think of commits as save points in a video game - you can always return to any previous commit if something goes wrong.

**Why commit?** Commits create a history of your work, making it easy to understand what changed, when, and why. They also make it possible to undo mistakes by reverting to earlier states.

### What is a Pull Request (PR)?

A **pull request** is a proposal to merge your changes from your fork back into the original repository. It's called a "pull" request because you're asking the maintainers to "pull" your changes into their project.

**Why pull requests?** PRs enable code review and discussion. Maintainers can review your code, suggest improvements, and ensure your changes align with the project's goals before merging them. This collaborative review process maintains code quality.

### What is Upstream vs Origin?

- **Origin**: The primary remote repository (for team members, this is the main PyBFS repo; for external contributors, this is your fork)
- **Upstream**: The original repository (only relevant for external contributors who have forked the project)

**For team members:** When you clone the main repository, `origin` refers to the PyBFS repository on GitHub. When you **pull** from origin, you're getting the latest changes. When you **push** to origin, you're uploading your branch to the main repository.

**For external contributors:** `origin` is your fork, and `upstream` is the main PyBFS repository. You **pull** from upstream to stay current and **push** to origin (your fork).

---

## Collaborative Workflow Overview

### For Core Team Members (Recommended)

If you have write access to the PyBFS repository:

1. **Clone the Repository** - Download the PyBFS repository to your local machine
2. **Create a Branch** - Make a new branch for your feature or fix
3. **Make Changes** - Implement your feature or bug fix
4. **Commit Changes** - Save your changes with descriptive commit messages
5. **Pull Latest Changes** - Sync with the main branch to avoid conflicts
6. **Push Your Branch** - Upload your branch to GitHub
7. **Submit a Pull Request** - Request to merge your changes into the main branch

### For External Contributors

If you don't have write access to the PyBFS repository:

1. **Fork the Repository** - Create your own copy of the project on GitHub
2. **Clone Your Fork** - Download your fork to your local machine
3. **Create a Branch** - Make a new branch for your feature or fix
4. **Make Changes** - Implement your feature or bug fix
5. **Commit Changes** - Save your changes with descriptive commit messages
6. **Pull Latest Changes** - Sync with the upstream repository to avoid conflicts
7. **Push to Your Fork** - Upload your changes to your fork on GitHub
8. **Submit a Pull Request** - Request to merge your changes into the main project

Both workflows keep the main repository clean while allowing you to work independently on features.

---

## How to Contribute: Choose Your Tool

Now that you understand the concepts and workflow, it's time to put them into practice. The sections below provide step-by-step instructions for each stage of the contribution process using three different development environments. Choose the environment that matches your preferred way of working:

- **Terminal (Bash/Command Line)** - For those who prefer command-line tools and want direct control over Git operations
- **PyCharm** - For Python developers who prefer JetBrains' integrated development environment
- **VS Code** - For developers who use Microsoft's popular, lightweight code editor

**Note:** Each section shows instructions for both team members (working directly with the repository) and external contributors (working with forks). Follow the instructions that match your access level.

---

## 1. Terminal (Bash/Command Line)

If you prefer to work directly from the command line, this section is for you. The terminal gives you precise control over Git operations and works on any operating system (macOS, Linux, or Windows with Git Bash). These commands form the foundation of all Git workflows, so understanding them will help you regardless of which tools you use later.

### Initial Setup

**For Team Members (with repository access):**

```bash
# Clone the main repository
git clone https://github.com/originalowner/pybfs.git
cd pybfs
```

**For External Contributors (without repository access):**

```bash
# First, fork the repository on GitHub (use the Fork button)

# Clone your fork
git clone https://github.com/yourusername/pybfs.git
cd pybfs

# Add upstream remote to sync with the main repository
git remote add upstream https://github.com/originalowner/pybfs.git
```

### Creating a New Branch

**For Team Members:**

```bash
# Make sure you're on main
git checkout main

# Pull latest changes from origin
git pull origin main

# Create and switch to a new branch
git checkout -b feature/your-feature-name
```

**For External Contributors:**

```bash
# Make sure you're on main
git checkout main

# Pull latest changes from upstream
git pull upstream main

# Create and switch to a new branch
git checkout -b feature/your-feature-name
```

### Making and Committing Changes

```bash
# Make your code changes, then check status
git status

# Add files to staging
git add .

# Commit with a descriptive message
git commit -m "Add your feature description"
```

### Pulling Latest Changes

**For Team Members:**

```bash
# Fetch and merge changes from origin main
git pull origin main

# Resolve any conflicts if they occur, then commit
```

**For External Contributors:**

```bash
# Fetch and merge changes from upstream main
git pull upstream main

# Resolve any conflicts if they occur, then commit
```

### Pushing Changes

```bash
# Push your branch to your fork
git push origin feature/your-feature-name

# If this is the first push, you may need:
git push -u origin feature/your-feature-name
```

### Submitting a Pull Request

**For Team Members:**

```bash
# Use GitHub CLI (if installed)
gh pr create --title "Your PR Title" --body "Description of changes"

# OR visit GitHub in your browser:
# Go to the repository on GitHub and click "Compare & pull request"
```

**For External Contributors:**

```bash
# Use GitHub CLI (if installed)
gh pr create --title "Your PR Title" --body "Description of changes"

# OR visit GitHub in your browser:
# Go to your fork on GitHub and click "Compare & pull request"
```

---

## 2. PyCharm

<img src="../images/pycharm.png" alt="PyCharm" width="80" style="float: left; margin-right: 15px; margin-bottom: 10px;">

[PyCharm](https://www.jetbrains.com/pycharm/) is a popular Integrated Development Environment (IDE) from JetBrains, designed specifically for Python development. It provides a rich graphical interface for Git operations, eliminating the need to memorize commands. PyCharm's built-in tools make it easy to visualize branches, compare changes, and manage pull requests without leaving your development environment. If you're already using PyCharm for coding, these steps will help you manage your contributions seamlessly.

### Initial Setup

**For Team Members (with repository access):**

1. **Clone the repository**:
   - Open PyCharm
   - Go to `Git` → `Clone...`
   - Enter the repository URL: `https://github.com/originalowner/pybfs.git`
   - Click `Clone`

**For External Contributors (without repository access):**

1. **Fork the repository** on GitHub (click Fork button)
2. **Clone your fork**:
   - Open PyCharm
   - Go to `Git` → `Clone...`
   - Enter your fork URL: `https://github.com/yourusername/pybfs.git`
   - Click `Clone`
3. **Add upstream remote**:
   - Go to `Git` → `Manage Remotes...`
   - Click `+` to add a new remote
   - Name: `upstream`
   - URL: `https://github.com/originalowner/pybfs.git`
   - Click `OK`

### Creating a New Branch

**For Team Members:**

1. **Switch to main branch**:
   - Click the branch name in the bottom-right corner
   - Select `main` from the list
2. **Pull latest changes**:
   - Go to `Git` → `Pull...`
   - Select `origin` as the remote
   - Select `main` as the branch
   - Click `Pull`
3. **Create new branch**:
   - Click the branch name in the bottom-right
   - Click `+ New Branch`
   - Name it: `feature/your-feature-name`
   - Make sure "Checkout branch" is checked
   - Click `Create`

**For External Contributors:**

1. **Switch to main branch**:
   - Click the branch name in the bottom-right corner
   - Select `main` from the list
2. **Pull latest changes**:
   - Go to `Git` → `Pull...`
   - Select `upstream` as the remote
   - Select `main` as the branch
   - Click `Pull`
3. **Create new branch**:
   - Click the branch name in the bottom-right
   - Click `+ New Branch`
   - Name it: `feature/your-feature-name`
   - Make sure "Checkout branch" is checked
   - Click `Create`

### Making and Committing Changes

1. **Make your code changes** in the editor
2. **View changes**:
   - Go to `Git` → `Commit` (or press `Cmd+K` / `Ctrl+K`)
3. **Stage and commit**:
   - Check the files you want to commit
   - Write a descriptive commit message
   - Click `Commit` (or `Commit and Push` to do both at once)

### Pulling Latest Changes

**For Team Members:**

1. **Pull from origin**:
   - Go to `Git` → `Pull...`
   - Select `origin` as the remote
   - Select `main` as the branch
   - Click `Pull`
2. **Resolve conflicts** if any appear in the Conflicts dialog
3. **Merge and commit** the changes

**For External Contributors:**

1. **Pull from upstream**:
   - Go to `Git` → `Pull...`
   - Select `upstream` as the remote
   - Select `main` as the branch
   - Click `Pull`
2. **Resolve conflicts** if any appear in the Conflicts dialog
3. **Merge and commit** the changes

### Pushing Changes

1. **Push to your fork**:
   - Go to `Git` → `Push...` (or press `Cmd+Shift+K` / `Ctrl+Shift+K`)
   - Verify the branch is pushing to `origin/feature/your-feature-name`
   - Click `Push`

### Submitting a Pull Request

**For Team Members:**

1. **Via PyCharm** (if GitHub plugin is enabled):
   - Go to `Git` → `GitHub` → `Create Pull Request`
   - Fill in the title and description
   - Click `Create Pull Request`
2. **Via Browser**:
   - PyCharm will show a notification with a link to create a PR
   - Click the link, or visit the repository on GitHub
   - Click `Compare & pull request`
   - Fill in details and click `Create pull request`

**For External Contributors:**

1. **Via PyCharm** (if GitHub plugin is enabled):
   - Go to `Git` → `GitHub` → `Create Pull Request`
   - Fill in the title and description
   - Click `Create Pull Request`
2. **Via Browser**:
   - PyCharm will show a notification with a link to create a PR
   - Click the link, or visit your fork on GitHub
   - Click `Compare & pull request`
   - Fill in details and click `Create pull request`

---

## 3. VS Code

<img src="../images/vscode.png" alt="VS Code" width="80" style="float: left; margin-right: 15px; margin-bottom: 10px;">

[Visual Studio Code (VS Code)](https://code.visualstudio.com/) is a lightweight, free, and highly extensible code editor from Microsoft that has become incredibly popular among developers. It offers a great balance between simplicity and power, with built-in Git support and a rich ecosystem of extensions. VS Code's Source Control panel provides an intuitive interface for Git operations while still giving you quick access to the terminal when needed. If you value a clean interface with excellent Git integration, VS Code is an excellent choice.

### Initial Setup

**For Team Members (with repository access):**

1. **Clone the repository**:
   - Open VS Code
   - Press `Cmd+Shift+P` / `Ctrl+Shift+P` to open Command Palette
   - Type and select `Git: Clone`
   - Enter the repository URL: `https://github.com/originalowner/pybfs.git`
   - Select a folder location
   - Click `Open` when prompted

**For External Contributors (without repository access):**

1. **Fork the repository** on GitHub (click Fork button)
2. **Clone your fork**:
   - Open VS Code
   - Press `Cmd+Shift+P` / `Ctrl+Shift+P` to open Command Palette
   - Type and select `Git: Clone`
   - Enter your fork URL: `https://github.com/yourusername/pybfs.git`
   - Select a folder location
   - Click `Open` when prompted
3. **Add upstream remote**:
   - Open the terminal in VS Code (`Ctrl+\`` or `View` → `Terminal`)
   - Run: `git remote add upstream https://github.com/originalowner/pybfs.git`

### Creating a New Branch

**For Team Members:**

1. **Switch to main branch**:
   - Click the branch name in the bottom-left corner
   - Select `main` from the dropdown
2. **Pull latest changes**:
   - Click the `...` menu in the Source Control panel
   - Select `Pull from...`
   - Choose `origin`
   - Select `main`
3. **Create new branch**:
   - Click the branch name in the bottom-left
   - Click `+ Create new branch...`
   - Name it: `feature/your-feature-name`
   - Press `Enter`

**For External Contributors:**

1. **Switch to main branch**:
   - Click the branch name in the bottom-left corner
   - Select `main` from the dropdown
2. **Pull latest changes**:
   - Click the `...` menu in the Source Control panel
   - Select `Pull from...`
   - Choose `upstream`
   - Select `main`
3. **Create new branch**:
   - Click the branch name in the bottom-left
   - Click `+ Create new branch...`
   - Name it: `feature/your-feature-name`
   - Press `Enter`

### Making and Committing Changes

1. **Make your code changes** in the editor
2. **View changes**:
   - Click the Source Control icon (branch icon) in the left sidebar
   - Or press `Ctrl+Shift+G`
3. **Stage and commit**:
   - Click `+` next to files to stage them (or `+` next to "Changes" to stage all)
   - Type a descriptive commit message in the message box
   - Click the `✓` checkmark (or press `Cmd+Enter` / `Ctrl+Enter`) to commit

### Pulling Latest Changes

**For Team Members:**

1. **Pull from origin**:
   - Click the `...` menu in the Source Control panel
   - Select `Pull from...`
   - Choose `origin`
   - Select `main`
2. **Resolve conflicts** if they appear:
   - VS Code will highlight conflicts in the editor
   - Click `Accept Current Change`, `Accept Incoming Change`, or `Accept Both Changes`
   - Stage and commit the merge

**For External Contributors:**

1. **Pull from upstream**:
   - Click the `...` menu in the Source Control panel
   - Select `Pull from...`
   - Choose `upstream`
   - Select `main`
2. **Resolve conflicts** if they appear:
   - VS Code will highlight conflicts in the editor
   - Click `Accept Current Change`, `Accept Incoming Change`, or `Accept Both Changes`
   - Stage and commit the merge

### Pushing Changes

1. **Push to your fork**:
   - Click the `...` menu in the Source Control panel
   - Select `Push`
   - Or click the sync icon (↻) in the bottom-left status bar
2. **Set upstream** if pushing for the first time:
   - VS Code will prompt you to publish the branch
   - Click `Publish Branch`

### Submitting a Pull Request

**For Team Members:**

1. **Via GitHub Pull Requests Extension** (if installed):
   - Install the "GitHub Pull Requests" extension if not already installed
   - Click the GitHub icon in the left sidebar
   - Click `Create Pull Request`
   - Fill in title and description
   - Click `Create`
2. **Via Browser**:
   - After pushing, VS Code may show a notification to create a PR
   - Click the notification link, or visit the repository on GitHub
   - Click `Compare & pull request`
   - Fill in details and click `Create pull request`

**For External Contributors:**

1. **Via GitHub Pull Requests Extension** (if installed):
   - Install the "GitHub Pull Requests" extension if not already installed
   - Click the GitHub icon in the left sidebar
   - Click `Create Pull Request`
   - Fill in title and description
   - Click `Create`
2. **Via Browser**:
   - After pushing, VS Code may show a notification to create a PR
   - Click the notification link, or visit your fork on GitHub
   - Click `Compare & pull request`
   - Fill in details and click `Create pull request`

---


## Questions?

If you have questions, please open an issue on GitHub.
