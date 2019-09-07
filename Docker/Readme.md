<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />  
  </head>
  <body>
<h1>DOCKER BASIC COMMANDS</h1>
<p>
  You can compile an image from this dockerfile with:<br>
<pre>docker build -t [image_name] .
eg. docker build -t rna_editing_protocol .</pre>
<br>
and run a container from it with:<br>
<pre>docker run -it [image_name] bash
eg. docker run -it rna_editing_protocol bash</pre>
    </p>
<p>
  <b>OR</b>
</p>
<p>
  download a pre-built image from dockerhub with:
  <pre>docker pull claudiologiudice/rna_editing_protocol:latest</pre>
  <br>
  and run a container from it with:
  <pre>docker run -it claudiologiudice/rna_editing_protocol:latest bash</pre>
  </p>
</body>
</html>
