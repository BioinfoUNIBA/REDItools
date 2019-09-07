You can compile your image from this dockerfile with:
docker build -t <image_name> .
eg. docker build -t rna_editing_protocol .

and run a container from it with:
docker run -it image_name bash
eg. docker run -it rna_editing_protocol bash


OR

download a pre-built image from dockerhub with:
docker pull claudiologiudice/rna_editing_protocol:latest

and run a container from it with:
docker run -it claudiologiudice/rna_editing_protocol:latest bash
