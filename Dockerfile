FROM python:2
USER root

MAINTAINER Guus Martijn Teunisse <gmteunisse@gmail.com>

RUN apt-get update \
	&& apt-get install -y \
	wget \
	glpk-utils \
	libglpk-dev \
	glpk-doc \
	&& pip install glpk \
	numpy \
	scipy \
	lxml \
	python-libsbml \
	swiglpk \
	&& rm -rf /var/lib/apt/lists/*

#RUN wget https://github.com/opencobra/cobrapy/archive/0.5.11.tar.gz \
#	&& tar -zxf 0.5.11.tar.gz \
#	&& cd cobrapy-0.5.11 \
#	&& pip install .
RUN pip install cobra==0.3.2

ADD /src /src
RUN chmod +x /src/assembly_model.py
ENV PATH="/src:${PATH}"
WORKDIR /src
ENTRYPOINT ["python", "/src/assembly_model.py"]
