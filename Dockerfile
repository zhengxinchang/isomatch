FROM rust:1.87-slim AS builder

WORKDIR /build
COPY . .
RUN cargo build --release

RUN mkdir /linux_build && cp target/release/isomatch /linux_build/

FROM scratch
COPY --from=builder /linux_build /linux_build
