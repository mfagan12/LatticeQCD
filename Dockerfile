FROM julia:latest
RUN apt update && apt-get update
RUN apt-get install -y python3 python3-pip \
	&& pip3 install numpy scipy jupyterlab
RUN julia -e 'using Pkg; Pkg.add("IJulia")'
