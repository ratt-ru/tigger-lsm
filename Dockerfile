FROM kernsuite/base:2

MAINTAINER gijsmolenaar@gmail.com

ADD . /tmp/tigger-lsm

RUN pip install /tmp/tigger-lsm

CMD /usr/local/bin/tigger-convert
