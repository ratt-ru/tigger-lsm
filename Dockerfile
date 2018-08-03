FROM kernsuite/base:4
RUN docker-apt-install \
    python-setuptools \
    python-numpy \
    python-scipy \
    python-astropy \
    python-astro-kittens \
    python-astlib \
    python-pip \
    python3-setuptools \
    python3-numpy \
    python3-scipy \
    python3-astropy \
    python3-astlib \
    python3-pip
RUN docker-apt-install git
ADD . /code
RUN pip3 install git+https://github.com/ska-sa/kittens.git@modernize
RUN pip install /code
RUN python2 /usr/local/bin/tigger-convert /code/test/3C147-HI6.refmodel.lsm.html /tmp/output.txt
RUN python2 /usr/local/bin/tigger-make-brick /code/test/3C147-HI6.refmodel.lsm.html /code/test/bla.fits
RUN python2 /usr/local/bin/tigger-tag /code/test/3C147-HI6.refmodel.lsm.html gijs
RUN python2 /usr/local/bin/tigger-restore -f /code/test/bla.fits /code/test/3C147-HI6.refmodel.lsm.html
RUN pip3 install /code
RUN python3 /usr/local/bin/tigger-convert -f /code/test/3C147-HI6.refmodel.lsm.html /tmp/output.txt
RUN python3 /usr/local/bin/tigger-make-brick /code/test/3C147-HI6.refmodel.lsm.html /code/test/bla.fits
RUN python3 /usr/local/bin/tigger-tag /code/test/3C147-HI6.refmodel.lsm.html gijs
RUN python3 /usr/local/bin/tigger-restore -f /code/test/bla.fits /code/test/3C147-HI6.refmodel.lsm.html
RUN echo "the next command should not print 1"
RUN wc -l /tmp/output.txt
