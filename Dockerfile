FROM kernsuite/base:3
RUN docker-apt-install python-pip
RUN docker-apt-install python-setuptools python-numpy python-scipy python-astropy python-kittens
RUN pip install astro-kittens astlib
ADD . /code
RUN pip install /code
RUN /usr/local/bin/tigger-convert /code/test/3C147-HI6.refmodel.lsm /tmp/output.txt
RUN /usr/local/bin/tigger-make-brick /code/test/3C147-HI6.refmodel.lsm /tmp/brick
RUN /usr/local/bin/tigger-tag /code/test/3C147-HI6.refmodel.lsm gijs
RUN /usr/local/bin/tigger-restore
RUN echo "the next command should not print 1"
RUN wc -l /tmp/output.txt
