Some help on SSH keys

Read the [chapter 4](https://learning.oreilly.com/library/view/bioinformatics-data-skills/9781449367480/ch04.html#chapter-04) from Vince Buffalo's book.

Also see information on setting up SSH keys via github.
[Generating SSH keys](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/)

```shell
$ ssh-keygen -t rsa -b 4096 -C "youremail"
# you will be prompted for a password - this is not your login password but a new password
# you can use for unlocking your keypair
```

Setting up ssh keys for passwordless login
==========================================
1. Generates a ssh key pair on your laptop
2. Need to copy ~/.ssh/id_rsa.pub to cluster.hpcc.ucr.edu

```shell
[YOURLAPTOP] $ scp  ~/.ssh/id_rsa.pub YOURNAME@cluster.hpcc.ucr.edu:.ssh/mylaptopkey.pub
[YOURLAPTOP] $ ssh YOURNAME@biocluster.ucr.edu
[hpcc] $ cd .ssh
[hpcc] $ cat mylaptopkey.pub >> authorized_keys
[hpcc] $ chmod 644 authorized_keys
[hpcc] $ logout
$ ssh YOURNAME@cluster.hpcc.ucr.edu
# should prompt you for your SSH key password
```

On your laptop enable local caching of the ssh key password by doing
```shell
$ ssh-add
password:
$ ssh  YOURNAME@biocluster.ucr.edu
# now no password requested!
```

Can edit ~/.ssh/config to setup aliases, preset the user name and simplify

```plain
$ cat ~/.ssh/config
ForwardX11 yes
ForwardX11Trusted yes
ForwardAgent yes

Host hpcc
 Hostname cluster.hpcc.ucr.edu
 User MYUSERNAME
 ServerAliveInterval 10
```

Can just use `ssh hpcc` now instead of `ssh MYUSERNAME@cluster.hpcc.ucr.edu`

Using Public SSH keys
====================
The public SSH keys from your laptop and one on biocluster can be
uploaded to github for easier checkin / checkout authentication via
SSH instead of HTTPS

Use this public key (`id_rsa.pub`) to your github account [https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account). Note that you can create multiple keys and have one pair for your laptop and another for HPCC. All of them can be added to github to enable check outs.
