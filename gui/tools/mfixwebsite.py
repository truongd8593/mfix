# -*- coding: utf-8 -*-
"""
This module provides utilities to log into the mfix website:
https://mfix.netl.doe.gov/

@author: Weberjm
"""
from __future__ import (division, absolute_import, print_function,
                        unicode_literals)

try:
    import cookielib
except:
    import http.cookiejar as cookielib
import urllib
try:
    # Python 2
    import urllib2
except:
    # Python 3
    import urllib.request as urllib2
import re
import sys
import os
import tarfile
import zipfile
import subprocess
from distutils.version import StrictVersion

class MfixWebsiteLoginError(Exception):
    """
    Base login error
    """
    pass

class MfixWebsiteFilenameError(Exception):
    """
    Base login error
    """
    pass

class MfixWebsite(object):
    """
    Class for loging into and downloading files from the mfix website.
    """
    def __init__(self, login=None, password=None, printprogress=True, path=None):
        """
        Initialize and login.
        """
        self.login = login
        self.password = password
        self.loggedin = False
        self.blocksize = 8192
        self.printprogress = printprogress
        self.path = path
        self.filelist = []
        self.mfixSourceLink = None
        self.mfixWindowsLink = None
        self.mfixWindowsVersion = None
        self.c3mLink = None
        self.c3mVersion = None
        self.guiLinkTar = None
        self.guiLinkZip = None

        # regex
        self.loginerror_regex = re.compile('<div id="login_error">(.*)<br />')
        self.fname_regex = re.compile('filename="(.*)"')
        self.dltable_regex = re.compile('alt="\[[ ]+\]"></td><td><a href="(.*)">(.*)</a>')
        self.windowsversion_regex = re.compile('[\d]+\-[\d]+')
        self.c3mversion_regex = re.compile('v([\d]+\.[\d]+\.[\d]+)')

        # build opener
        self.cj = cookielib.CookieJar()
        self.cj.clear()
        self.opener = urllib2.build_opener(
            urllib2.HTTPRedirectHandler(),
            urllib2.HTTPHandler(debuglevel=0),
            urllib2.HTTPSHandler(debuglevel=0),
            urllib2.HTTPCookieProcessor(self.cj)
        )
        self.opener.addheaders = [
            ('User-agent', ('mfixGui/1.0')),
        ]

        # get cookies
        self._fillCookieJar()

        # make directory if it doesn't exisit
        if self.path and not os.path.exists(self.path):
            os.mkdir(self.path)

        # login
        if self.login and self.password:
            self.loginWebsite()

    def setdownloadpath(self, path):
        self.path = path

        if self.path and not os.path.exists(self.path):
            os.mkdir(self.path)

    def _fillCookieJar(self):
        """
        Fill the cookie jar.
        """
        self.opener.open("https://mfix.netl.doe.gov/wp-login.php")

    def loginWebsite(self):
        """
        Handle login.
        """
        login_data = urllib.urlencode({
            'log' : self.login,
            'pwd' : self.password,
            'rememberme' : 'forever',
            'wp-submit' : 'Log In',
            'testcookie' : '1',
            'redirect_to':'https://mfix.netl.doe.gov/wp-admin/',
        })

        response = self.opener.open("https://mfix.netl.doe.gov/wp-login.php", login_data)

        for line in response.readlines():
            logErr = self.loginerror_regex.findall(line)
            if logErr:
                logErr = logErr[0].strip()
                break

        if logErr:
            logErr = re.sub('<.+?>','', logErr)
            raise MfixWebsiteLoginError(logErr)
        else:
            self.loggedin = True

        self.crawlDownload()

    def crawlDownload(self):
        self.mfixLinks = self.extractLinksFromtable('https://mfix.netl.doe.gov/download/mfix/')
        self.c3mLinks = self.extractLinksFromtable('https://mfix.netl.doe.gov/download/c3m/')
        self.mfixGuiLinks = self.extractLinksFromtable('https://mfix.netl.doe.gov/download/mfix/mfix-gui/')

        # find mfix links
        for mfixlink in self.mfixLinks:
            if 'windows' in mfixlink[1]:
                self.mfixWindowsLink = mfixlink[0]

                if self.windowsversion_regex.findall(mfixlink[1]):
                    self.mfixWindowsVersion = self.windowsversion_regex.findall(mfixlink[1])[0]
                elif self.windowsversion_regex.findall(mfixlink[0]):
                    self.mfixWindowsVersion = self.windowsversion_regex.findall(mfixlink[0])[0]
                else:
                    self.mfixWindowsVersion = None
            elif mfixlink[1] == 'mfix.tar.gz':
                self.mfixSourceLink = mfixlink[0]

        # c3m links
        for c3mlink in self.c3mLinks:
            version = self.c3mversion_regex.findall(c3mlink[1])
            if version:
                version = version[0]
                if self.c3mVersion:
                    if StrictVersion(self.c3mVersion) < StrictVersion(version):
                        self.c3mVersion = version
                        self.c3mLink = c3mlink[0]
                else:
                    self.c3mVersion = version
                    self.c3mLink = c3mlink[0]

        # mfix gui links
        for guiLink in self.mfixGuiLinks:
            if guiLink[1].endswith('.tar.gz'):
                self.guiLinkTar = guiLink[0]
            elif guiLink[1].endswith('.zip'):
                self.guiLinkZip = guiLink[0]

    def extractLinksFromtable(self, url):

        response = self.opener.open(url)

        linkList = []
        for line in response.readlines():
            fileLink = self.dltable_regex.findall(line)
            if fileLink:
                fileLink = fileLink[0]
                linkList.append(['/'.join([url,fileLink[0]]), fileLink[1]])

        return linkList

    def download(self, url, fname=None, name=None):
        if not self.loggedin:
            raise MfixWebsiteLoginError('Not logged into https://mfix.netl.doe.gov/')

        self.response = self.opener.open(url)

        meta = self.response.info()
        meta_func = meta.getheaders if hasattr(meta, 'getheaders') else meta.get_all
        meta_length = meta_func("Content-Length")
        self.file_size = None
        if meta_length:
            self.file_size = int(meta_length[0])

        self.file_size_dl = 0
        self.downloadprogress = 0.0

        if self.printprogress:
            if fname is None:
                if 'content-disposition' in meta:
                    fname = self.fname_regex.findall(meta['content-disposition'])[0]
                else:
                    raise MfixWebsiteFilenameError('Could not find a filename to use.')

            if self.path is not None:
                fname = os.path.join(self.path, fname)

            self.filelist.append(fname)
            self.fname = fname

            with open(fname, 'wb') as f:
                while True:
                    chunk = self._downloadchunk()
                    if chunk is not None:
                        f.write(chunk)
                    else:
                        break


                    p = int(self.downloadprogress/5)
                    status = 'Downloading'
                    if name:
                        status+=' '+name
                    status += ': ['+'#'*p+' '*(20-p)+']'
                    print(status+'{0:6.0f}%'.format(self.downloadprogress), end='\r')
                print()


    def _downloadchunk(self):
        '''
        Download a chunk of the file. chunk size: self.blocksize
        '''
        buffer = self.response.read(self.blocksize)
        if not buffer:
            return None

        self.file_size_dl += len(buffer)
        if self.file_size:
            self.downloadprogress = self.file_size_dl/self.file_size*100

        return buffer

    def downloadMfixReleaseTarball(self,):
        '''
        Download mfix tarball (.tar.gz).
        '''
        if self.mfixSourceLink:
            self.download(self.mfixSourceLink, fname='mfix.tar.gz', name='mfixtarball')

    def downloadMfixReleaseWindows(self,):
        '''
        Download windows binary.
        '''
        if self.mfixWindowsLink:
            self.download(self.mfixWindowsLink, fname='mfix.zip', name='mfixwindows')

    def downloadC3M(self,):
        '''
        Download C3M installer.
        '''
        if self.c3mLink:
            self.download(self.c3mLink, fname='c3minstaller.exe', name='c3m')

    def downloadMfixGui(self, comp=None):
        '''
        Download the mfix gui.
        '''

        if comp is None:
            if os.name == 'nt':
                comp = 'zip'

        if comp =='zip':
            if self.guiLinkZip:
                self.download(self.guiLinkZip, fname='mfixgui.zip', name='mfixgui')
        else:
            if self.guiLinkTar:
                self.download(self.guiLinkTar, fname='mfixgui.tar.gz', name='mfixgui')


    def extract(self, fname, extractpath='.'):

        if fname.endswith('.tar.gz'):
            tar = tarfile.open(fname, 'r:gz')
            tar.extractall(extractpath)
            tar.close()

            rename = os.path.basename(fname)
            for i in range(3):
                rename = os.path.splitext(rename)[0]

            while os.path.exists(rename):
                rename += '_1'

            os.rename(os.path.join(extractpath, 'mfix'),
                      os.path.join(extractpath, rename),
                      )

        elif fname.endswith('.zip'):
            zipf = zipfile.ZipFile(fname)
            extractpath = os.path.basename(os.path.splitext(fname)[0])
            extractpath = os.path.join(os.path.dirname(fname), extractpath)

            zipf.extractall(extractpath)
            zipf.close()

    def install(self, fname):
        subprocess.Popen(fname)

if __name__ == '__main__':
    import getpass

    if '-help' in sys.argv:
        print('Arguments:',
              '',
              'Download Options',
              '================',
              '-downloadall            download the mfix suite',
              '-downloadwindows        download the mfix suite for windows',
              '-downloadlinux          download the mfix suite for linux',
              '-mfixtarball            download tarball',
              '-mfixwindows            download windows binary',
              '-mfixgui                download mfix gui; .zip on windows, .tar.gz on unix',
              '-c3m                    download c3m installer',
              '-downloadpath <path>    <path> to save files too, default is cwd',
              '',
              'Tools',
              '=====',
              '-extract                extract compressed files',
              '-install                install *.exe on windows',
              '-help                   print this help',
              sep='\n')
        sys.exit()

    usr = raw_input('User Name: ')
    pwd = getpass.getpass()

    try:
        web = MfixWebsite(usr, pwd)
    except MfixWebsiteLoginError as e:
        print(e)
        sys.exit()

    if '-downloadpath' in sys.argv:
        path = sys.argv[sys.argv.index('-downloadpath')+1]

        if not os.path.exists(path):
            print('Path does not exsist: {}'.format(path))
            sys.exit()
    else:
        path = './'

    web.setdownloadpath(path)

    try:
        if '-downloadall' in sys.argv:
            web.downloadMfixReleaseTarball()
            web.downloadMfixReleaseWindows()
            web.downloadMfixGui()
            web.downloadC3M()

        elif '-downloadwindows' in sys.argv:
            web.downloadMfixReleaseWindows()
            web.downloadMfixGui()
            web.downloadC3M()

        elif '-downloadlinux' in sys.argv:
            web.downloadMfixReleaseTarball()
            web.downloadMfixGui()

        else:
            if '-mfixtarball' in sys.argv:
                web.downloadMfixReleaseTarball()

            if '-mfixwindows' in sys.argv:
                web.downloadMfixReleaseWindows()

            if '-mfixgui' in sys.argv:
                web.downloadMfixGui()

            if '-c3m' in sys.argv:
                web.downloadC3M()

    except urllib2.HTTPError as e:
        print(e)
        print('Please try again in a few minutes.')
        sys.exit()

    if '-extract' in sys.argv:
        for fname in web.filelist:
            if fname.endswith('.zip') or fname.endswith('.tar.gz'):
                print('Extracting {} ...'.format(os.path.basename(fname)), end='')
                web.extract(fname, extractpath=path)
                print('Done!')

    if '-install' in sys.argv:
        if os.name == 'nt':
            for fname in web.filelist:
                if fname.endswith('.exe'):
                    web.install(fname)
        else:
            print('Can only install on windows.')
