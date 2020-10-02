Git Tutorial
===========
1. Create account on github.com - see [https://help.github.com/](https://help.github.com/)
    * Also see [GitHub book](http://proquest.safaribooksonline.com/9781449367480/ch05_shared_central_github_html)

2. Login

3. Create new repository 'mytestrepo' or click on homework link to automatically create your homework repository for the specific assignment

4. Add a README.md file. This ends in 'md' because we usually write
these in [markdown](https://guides.github.com/features/mastering-markdown/)
which is a structured text format that can be rendered in nice ways in
HTML and PDF (like this guide). You can also make it a plain text file
and use '.txt'.

5. Upload files to the repository online.
6. Checkout the code.
```shell
$ git clone https://github.com/YOURID/mytestrepo.git
```
7. Edit / Update a file
```shell
$ nano solution.sh
# or some other tool (e.g. emacs, textedit, vi, ...)
```

8. Commit these changes to your LOCAL repository
```shell
$ git commit -m "I updated the solution"
```

9. Save this changed repository in the github server. (eg 'the
cloud'). This lets others see the file changes and sync their
repository with yours if you are doing a collaborative project or for
the instructor to see your submitted solutions. Also saves a copy, you
or others can then checkout a clone of this repository on other
computers or even just different folders if you need to
```shell
$ git push
# git push will request your github ID and github password
```
* See [SSH_keys](SSH_keys.html) info about how to make SSH keys and [Github info on adding SSH keys](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/) to your account so you don't have to enter password each time.
10. Add another new file (e.g. another script)
```shell
$ nano newfile.py
$ git add newfile.py
$ git commit -m "I added this new file" newfile.py
$ git push
```
