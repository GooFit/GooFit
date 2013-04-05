README for GooFit developers
============================

Docs
----

Local
+++++

We use `Doxygen <http://www.doxygen.org>`_ to generate GooFit API docs.
To generate html documentation locally on your machine and view it run
``doxygen`` in the top-level folder of the GooFit repo and then open
``html/index.html`` in your browser.

Push Updated Version Online
+++++++++++++++++++++++++++

We use `Github Pages <http://pages.github.com>`_ to serve GooFit API docs
online at `http://GooFit.github.com/GooFit/`_. If you want to update
the online version with your local ``html`` folder contents for the first
time, execute these steps::

	rm -r html # remove html folder if it existed
	mkdir html
	cd html
	git clone git@github.com:GooFit/GooFit.git .
	git checkout gh-pages

Now when you are in the ``html`` folder you are in a separate clone of the
``GooFit`` repo in the ``gh-pages`` branch. If you ``cd ..`` into the top-level
folder you are back in the original clone of the ``GooFit`` repo where
you do everything else but push documentation to ``gh-pages``.

Now that you are set up, it's easy to update the online docs::

	cd <top-level GooFit folder>
	doxygen # generate docs in html folder
	cd html
	# If you want to review the changes
	git status
	git diff
	git add .
	git commit -m 'update docs'
	git push origin gh-pages

Then go to `https://github.com/GooFit/GooFit/tree/gh-pages`_ or
`http://GooFit.github.com/GooFit/`_ to review your changes.
