# Contributing
We appreciate your consideration in supporting the development of Optimus Bind. The repository is well underway with a defined program outline. If you wish to contribute now is a great time to get involved!

## State of Development:
The current state of development and full listing of tasks can always be viewed at [Trello](https://trello.com/invite/b/V94BBx1d/4550ff50fe61eb27b8d304da57b00fe8/optimus-bind). Please coordinate with us on slack if you wish to contribute.

<strong>Milestones Underway:</strong>
<details>
<summary>Scoring:</summary>

-  IAlign.pl and ISAlign.pl first drafts are being finished up.
-  work on PSI-Blast and MSA is in progress!
    <img src="https://lh3.googleusercontent.com/Yl9t4Qm6OU5cUxEEDBkifFS3pSrJe73Kk_oLAfhPB99GxBd1KaSOSjo9c-bv-O-8kmi5DYyCCNoKGXh5NGaRfCzq6JrSK3YhfMHpI3KLYO-Lcd-AvNlTN1H-o7lxBqGa5wQ_BimykjymOLeUygJWYmew1OJ0AuprrLyfdazKZuhmYwKn2Z9JC-7f26N-LNYgaS-5Oa54J9C8N25dQWj7-kXjF5Cb0XwmgWJr4ceMlv8K2RH_3IIfaIiQSG2aj4ViigcHLrN7hd9BT4ZtMW5fWKHqTK1kinWa2YuSnJcmys6urgWEmjxb1XMo5YkBdquvKkYyUjTE-ckgHNh4uGhk4gjXWAHxWgo9Lir0KzRMRuzQw5Cy6DcfDf8aIosOebm0bFdFrK-3M4bETHZ1EPKKEh1fEPPpEVSpa2znHf9cpY0zEm7XjG7o-tVXw1GREqYg8HUBjNG3PY2fZaIP2QCuRb1e--UC9KVK13CosHqH4B5yi1uaiPyIbRDHjBO3XoCrEB7brNvNJTGvZpgiSgZ3Gy3_h2N8Xojk4tTwYmsHW4YydfDt6jNLhgBzWmmSJuCWyfYU_xi-J_78xG2LZpijTPlsPKGcG47bGnNant9ArlUJiA008DT9OH1pknzrGzxQYCp2QBJVa-IkeBH-fL2b9cnWudDaBbsFqQQ0T2I0-myFcwx1DR2S5HAMDDJ0GKhmqk1_5nBbop6xjyYj5uRIxNF0=w1420-h966-no"> 
</details>

<details>
<summary> Machine learning:</summary>

- Scoring required first…
</details>

<details>
<summary> Holistic Integration:</summary>

-   Testing is in its infancy
-   Check 3rd party requirements
-   Code cleanup (just revising and reviewing)
</details>

## Getting Started:
As an open source project, we encourage a diverse set of backgrounds. As such, it is anticipated that you may not be familiar with all the recourses we utilize. If you require help getting set up, please do not hesitate to ask!

###  1) Accounts:
- *Repository:*  [GitHub](https://github.com/tcardlab/optimus_bind_sample)
- *To-do list:*  [Trello](https://trello.com/invite/b/V94BBx1d/4550ff50fe61eb27b8d304da57b00fe8/optimus-bind) (Click to join)
- *Chat room:*  [Slack](https://bioscienceclub.slack.com/messages/CHK7D10MN/details/) (Upon request)

### 2) Environment and Installs:
- Refer to [README.md](https://github.com/tcardlab/optimus_bind_sample#installing)

### 3) Workflow (GitFlow): 
<details open><summary>Toggle Details</summary>

If you are new to [GitFlow](https://www.youtube.com/watch?v=aJnFGMclhU8&t=194), just think of it as a naming convention for branches. Doing so helps to simultaneously manage the development, release, and maintenance of the repo. 

Video explaining this [Workflow in practice](https://www.youtube.com/watch?v=Lj_jAFwofLs) (~5min on 1.5 speed).

<details>
<summary>Branching Conventions:</summary>

*Core branches:* `master & develop` (do not commit to or duplicate these)

*Work branches:* `branchType/your-branch-title` (branchTypes blow)

<u>Master branches<u>
- `develop`
- `hotfix/... `

<u>Develop branches:<u>
 - `feature/... `
 - `bugfix/... `
-  `release/...`
<br></details>

<details>
<summary>Reasoning:</summary>

- Explicit/extensive workflow yields consistent practice & organized contributions  
- Highly segmented development lends itself to well organized feature branches
- Not dependent upon the rapid development of continuous integration
</details>
</details>

### 4) On using git:
<details>
<summary><i>Git-Beginner (GUI):</i></summary>

If you have no Git/GitFlow experience and just want to get to work:
- [Github Desktop](https://desktop.github.com/)
	- [Setup](https://help.github.com/en/desktop/getting-started-with-github-desktop/setting-up-github-desktop)
	- The program itself is inherently directive, very easy to pick up.
- [SourceTree](https://www.sourcetreeapp.com/)
	- [GitFlow-Intro](https://medium.com/@budioktaviyans/how-to-make-a-git-flow-using-sourcetree-20ab77fe6813)

As a beginner, there is no need to concern yourself with the specifics of GitFlow, just the proper naming conventions.

Fork 
https://github.com/tcardlab/optimus_bind_sample/fork

<br></details>

<details>  
<summary><i>Git-Veteran (CLI):</i></summary>

If you are new to command line git, review [here](https://dont-be-afraid-to-commit.readthedocs.io/en/latest/git/commandlinegit.html).
<h4>Feature Contribution:</h4>

Join as a contributor or [fork](https://github.com/tcardlab/optimus_bind_sample/fork).
```
0) Clone appropriate repo and enter directory

1) Open [develop] branch
    $ git checkout develop  
    
2) Create new branch off [develop] (-b)
   Use one branch per feature / fix
    $ git checkout -b feature/your-branch-name

3) Commit changes in relevant chunks as work proceeds 
    $ git commit -am 'short commit description'

4) Share code online
    $ git push origin your-branch-name

5) Submit changes (Performed on GitHub)
  # If Contributor:
    Submit pull request to [develop] branch
	
  # If Forked:
    Submit pull request of [your-branch-name] 
    to original repo at [develop]
```
[NOTE: Commits may occur after pull request](https://help.github.com/en/articles/committing-changes-to-a-pull-request-branch-created-from-a-fork). After code review, if changes need to be made, the pull request will automatically update with new commits. If you are to clone a forked directory, make sure you clone to a different directory than where you keep the original clone.

[GitFlow extensions](https://danielkummer.github.io/git-flow-cheatsheet/) may provide useful for anything more complex. Primarily in regard to more complex branch management and releases. 
(Extensions currently maintained at: [GitFlow-avh](https://github.com/petervanderdoes/gitflow-avh))
<br></details>
[Pull request checklist](https://github.com/kylelobo/The-Documentation-Compendium/blob/master/en/PULL_REQUEST_TEMPLATE.md)

### 5) Style and Structure:
<details>
<summary><i>Codebase Structure:</i></summary>

Structure and organization is provided by [CookieCutter DataSci](https://drivendata.github.io/cookiecutter-data-science).

See [ReadMe](https://github.com/tcardlab/optimus_bind_sample#project-organization) for project outline.
</details>

<details>
<summary><i>Python:</i></summary>

Python submissions should follow [PEP8 guidelines](https://www.python.org/dev/peps/pep-0008/). Use a [PEP8 checker](http://pep8online.com/) or sufficient linter to ensure validity.

*Exception:* comment spacing may be multiple of 2 as it is useful for collapsing section within functions in IDE's.
 <br></details>

<details>
<summary><i>Commit Messages:</i></summary>

[Review  best practice here](https://gist.github.com/robertpainsi/b632364184e70900af4ab688decf6f53#file-commit-message-guidelines-md)

Please make use of commit [descriptions](https://stackoverflow.com/questions/9562304/git-github-commit-with-extended-message-description):
-   Describe why a change is being made.
-   How does it address the issue?
-   Do not assume the reviewer understands what the original problem was.
-   Do not assume the code is self-evident/self-documenting.
</details>

### 6) Testing:
<details>
<summary></summary>
IDK... refer to read me.~
lol
<br></details>


<br><br>
## While under development please refer the following topics to slack discussion:
### Testing
### File bug report
### Suggest new features
