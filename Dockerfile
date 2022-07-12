FROM bitnami/minideb:bullseye AS build

ADD https://github.com/andersgs/ph-cue/archive/refs/tags/v0.2.0.tar.gz /tmp/ph-cue.tar.gz

RUN apt update -y && apt install -y curl build-essential

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y 

ENV PATH $PATH:/root/.cargo/bin

RUN cd /tmp && tar -xzvf ph-cue.tar.gz && cd ph-cue-0.2.0 && cargo build -r

FROM bitnami/minideb:bullseye AS run

COPY --from=build /tmp/ph-cue-0.2.0/target/release/ph-cue /usr/local/bin/ph-cue

