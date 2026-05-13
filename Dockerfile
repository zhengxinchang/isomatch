FROM rust:1.92-slim AS builder

WORKDIR /build
COPY . .
RUN cargo build --release

RUN mkdir /linux_build && cp target/release/isomatch /linux_build/
