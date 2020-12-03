# Genome Browsers for interacting with Genomic data

Many different browser environments

1. [NCBI Genome](https://www.ncbi.nlm.nih.gov/genome) [Maize](https://www.ncbi.nlm.nih.gov/genome/12) [Maize chromosome](https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_902167145.1)
2. [Ensembl](https://ensembl.org)
3. [FungiDB](https://fungidb.org)
4. [Mycocosm](https://mycocosm.jgi.doe.gov/mycocosm/home)
4. [Wormbase](https://wormbase.org/)
5. [FlyBase](https://flybase.org/)
6. [Saccharomyces Genome Database](https://yeastgenome.org/)
7. [TAIR](http://arabidopsis.org)
8. [Gramene](https://www.gramene.org/) - Plant Comparative Resources
9. [Phytozome](https://phytozome.jgi.doe.gov/) - Plant Comparative Genomics portal
5. [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)

Setting up your own - JBrowse - Genome Browser
=====

(some of this was first written by: Tania Kurbessoian [@tania-k](https://github.com/tania-k) )

To visualize genome annotation combined with Epigenomic, Transcriptomic, or Variant data you you want to visualize them onto a genome browser.  [JBrowse](https://jbrowse.org) provides an easy to setup tool for this visuzalition.

There is substantial [documentation](https://jbrowse.org/docs/installation.html) that describes installation and add-on features. The [FAQ](https://jbrowse.org/docs/faq.html) is also incredibly helpful.

# Setup JBrowse on HPCC

These steps will show you how to setup JBrowse on HPCC with some already installed systems to make it easier for you

### Configure your HPCC account to be able to share via HTTP / Web

First you need to configure your account to be able to share data via the web.

Follow the [directions on the HPCC manual](http://hpcc.ucr.edu/manuals_linux-cluster_sharing.html#sharing-files-on-the-web) so that you can configure your home folder `~/.html` to be able to serve up data. If you do not want to make everything in this folder public you can use some simple strategies to enable a password protected space by [creating a `.htaccess`](http://hpcc.ucr.edu/manuals_linux-cluster_sharing.html#password-protect-web-pages) file or you can use a web form of security - security through obscurity.

```
mkdir ~/.html/private
chmod 711 ~/.html/private
mkdir ~/.html/private/MYSECRETPROJECT
# install jbrowse or data or others in ~/.html/private/MYSEKR1TPROJECT
```

On the web a user browsing will not have permission to see `https://cluster.hpcc.ucr.edu/~YOURUSERNAME/private` but the URL `https://cluster.hpcc.ucr.edu/~YOURUSERNAME/private/MYSEKR1TPROJECT` will be visible.  This is called security by obscurity, if you generate a long random string instead of `MYSEKR1TPROJECT` it would be hard to guess it (though note this not really secure since anyone reading network traffic and see the string would now know the the URL to go to).

Generally if you want to protect the data, setup a `.htaccess/.htpasswd` to require logging in. [Making an htpasswd file](https://hpcc.ucr.edu/manuals_linux-cluster_sharing.html#password-protect-web-pages)

## Setting up your own copy of JBrowse software

The next directions are specific to the UCR HPCC. These instructions use an already build conda environment which you can link to:
```
# startup an interactive job on the cluster
srun -p short -N 1 -n 4 --mem 16gb --pty bash -l
mkdir -p ~/bigdata/jbrowse2
cd ~/.html/private
ln -s ~/bigdata/jbrowse2 .
cd ~/.html/private

module load jbrowse/2
jbrowse create SARS-CoV-2
```

If you are going to support multiple JBrowse environments you only need to have a custom data folder. So you can symlink to all the files within the jbrowse checkout and then make a separate data folder too. Otherwise you need to make sure you have a separate custom jbrowse checkout for each project you are supporting.
