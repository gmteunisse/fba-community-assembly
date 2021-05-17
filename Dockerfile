FROM python:2
USER root

MAINTAINER Guus Martijn Teunisse <gmteunisse@gmail.com>

RUN apt-get update \
	&& apt-get install -y \
	glpk-utils \
	libglpk-dev \
	&& pip install \
	numpy \
	scipy \
	lxml \
	python-libsbml \
	&& rm -rf /var/lib/apt/lists/*

RUN pip install cobra==0.3.2

ADD /src /src
RUN chmod +x /src/assembly_model.py
ENV PATH="/src:${PATH}"
WORKDIR /src
ENTRYPOINT ["python", "/src/assembly_model.py"]
