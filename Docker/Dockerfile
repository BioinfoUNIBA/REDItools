#Download base image centos latest 
FROM centos

#Dockerfile Mantainer
LABEL mantainer="clalogiudice@gmail.com"

#Update the centos software with yum package-manager
RUN yum update -y && yum clean all

#Install git, wget and nano package
RUN yum -y install git wget nano && yum clean all

#Clone Nature_protocol Git repository
RUN git clone https://github.com/BioinfoUNIBA/REDItools

WORKDIR "/REDItools/NPscripts/" 

#Install miniconda with conda packages required by the nature_protocol
RUN chmod +x conda_pckg_installer_docker.py
RUN ./conda_pckg_installer_docker.py
ENV PATH /miniconda2/bin:$PATH
RUN echo "source activate nature_protocol" >> ~/.bashrc

#PREPARE NATURE_PROTOCOL environment
WORKDIR "/"
RUN echo "python ./REDItools/NPscripts/download-prepare-data-NP_docker.py" >> /root/.bashrc
