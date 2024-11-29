FROM harbor.ngs-ai.com/public/debian:bullseye-slim
RUN apt-get update && apt-get install -y procps
COPY target/x86_64-unknown-linux-gnu/release/orfan /usr/bin/
RUN /usr/bin/orfan --help
