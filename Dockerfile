FROM oraclelinux:7 AS builder

RUN yum -y update && yum clean all && \
    yum -y install \
        wget curl gcc gcc-c++ make \
        zlib-devel \
        which git

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo
ENV PATH=/usr/local/cargo/bin:$PATH

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y && \
    rustup default stable

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release

FROM amazonlinux:2
WORKDIR /linux_build
COPY --from=builder /build/target/release/isomatch ./isomatch
ENTRYPOINT ["./isomatch"]
