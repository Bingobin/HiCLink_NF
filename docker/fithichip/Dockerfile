FROM aylab/fithichip@sha256:bbb2cbcfc4bc61123d5100b111b2a648100da5fe8d24236ca0e4a3b7d8162fb4

# Keep the official runtime, but supply the current official scripts and missing dependencies.
COPY Miniconda3-py310_25.1.1-2-Linux-x86_64.sh /tmp/python-installer.sh
RUN bash /tmp/python-installer.sh -b -p /opt/fithichip-python
ENV PATH="/opt/fithichip-python/bin:${PATH}"
RUN python -m pip install --no-cache-dir numpy==1.26.4 pandas==1.5.3 scipy==1.15.3 h5py==3.11.0 cooler==0.10.4 networkx==2.8.8
COPY data.table_1.14.8.tar.gz /tmp/data.table.tar.gz
RUN R CMD INSTALL /tmp/data.table.tar.gz
COPY FitHiChIP/ /FitHiChIP/
COPY bedToBigBed /usr/local/bin/bedToBigBed
COPY interact.as /FitHiChIP/interact.as
RUN chmod +x /usr/local/bin/bedToBigBed
LABEL org.opencontainers.image.source="https://github.com/ay-lab/FitHiChIP" \
      org.opencontainers.image.revision="0ea1ac21be870908c672316ffbb630189dc6fae2"
ENTRYPOINT ["/bin/bash"]
