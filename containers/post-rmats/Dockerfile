FROM continuumio/miniconda:4.6.14

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

RUN pip install pandas

ADD ./sampleCountsSave.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/sampleCountsSave.sh
ADD ./create-matrix.py /usr/local/bin/
RUN chmod +x /usr/local/bin/create-matrix.py
ADD ./create_matrices_from_files.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/create_matrices_from_files.sh
ADD ./normalize_matrices_from_files.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/normalize_matrices_from_files.sh
