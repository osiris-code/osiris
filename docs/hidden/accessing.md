# Accessing the code

The Osiris open-source code is hosted on github at [https://github.com/osiris-code/osiris].

## Downloading and maintaining the code (NEEDS UPDATE)

Here's a short primer on working with the OSIRIS subversion repository:

### 1. Getting a fresh Osiris copy

To get a fresh OSIRIS copy from the development branch (3.0) to the current directory (this is called 'checking out' the code):

```text
 svn co https://osiris.ist.utl.pt/svn/branches/dev_3.0 .
```

Note the dot at the end of the command. When prompted just give the login credentials you were issued.

### 2. Update your existing copy

To update your existing copy with the latest changes do the following in the directory where you placed the code:

```text
svn update
```

And that's it! More information about using subversion can be found in the Subversion book available (for free) online at <http://svnbook.org/> and/or by issuing the 'svn help' command.

Any changes to the Osiris source code will be posted automatically to the osiris-users mailing list, so remember to keep your version up to date. This works very well even in case you made your own changes to the code: when updating, subversion will merge the changes for you. Also, remember that it is always possible (and easy) to revert to a previous version if something got broken.
