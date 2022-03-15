FROM alpine:3.14
COPY root/* /root/
WORKDIR /root
RUN apk update
RUN apk add build-base
RUN apk add openssh
RUN apk add --upgrade openmpi openmpi-dev
CMD [ "sh", "-c", "./script.sh" ]
