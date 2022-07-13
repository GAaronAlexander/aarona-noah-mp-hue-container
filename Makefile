# https://www.gnu.org/software/make/manual/html_node/Makefile-Conventions.html

SHELL = /bin/bash

.ONESHELL:
.SUFFIXES:

AWS_REGION = us-east-1
AWS_REGISTRY = 660097632732.dkr.ecr.$(AWS_REGION).amazonaws.com
GIT_REGISTRY = registry.gitlab.com
ORG = jupiterintel
REPO = aarona-noah-mp-hue
AWS_PROFILE ?= default
VERSION ?= latest
ARGS ?=
BUILD_ARGS ?=
IMAGE = $(ORG)/$(REPO):$(VERSION)


build:
	@git rev-parse HEAD > version && \
	DOCKER_BUILDKIT=1 && export DOCKER_BUILDKIT && \
	  docker build -f Dockerfile -t $(IMAGE) $(BUILD_ARGS) --ssh=default . && \
	rm version

clean:
	@rm -rf report.xml
	@find . -type d -name '__pycache__' -exec rm -rf {} +
	@find . -type d -name '*pytest_cache*' -exec rm -rf {} +
	@find . -type f -name "*.py[co]" -exec rm -rf {} +

push: build
	@docker tag $(IMAGE) $(GIT_REGISTRY)/$(IMAGE)
	@docker push $(GIT_REGISTRY)/$(IMAGE)

release: build
	@aws ecr get-login-password --region $(AWS_REGION) | docker login --username AWS --password-stdin 660097632732.dkr.ecr.$(AWS_REGION).amazonaws.com
	@docker tag $(IMAGE) $(AWS_REGISTRY)/$(IMAGE)
	@docker push $(AWS_REGISTRY)/$(IMAGE)

run:
	@AWS_ACCESS_KEY_ID=$$(aws --profile $(AWS_PROFILE) configure get aws_access_key_id)
	@AWS_SECRET_ACCESS_KEY=$$(aws --profile $(AWS_PROFILE) configure get aws_secret_access_key)
	@docker run --rm -it \
		-e AWS_ACCESS_KEY_ID="${AWS_ACCESS_KEY_ID}" \
		-e AWS_SECRET_ACCESS_KEY="${AWS_SECRET_ACCESS_KEY}" \
		$(IMAGE) $(ARGS)

test: build
	@python /home/jupiter/model/noahmp/generate-era5-boundary-conditions.py --start-date 2018-04-01 --end-date 2018-10-01 --freq 1H --save-location s3://jupiter-intern-projects/aarona/noah-mp-hue/milwaukee/2018-04-01_2018-10-01/ICBC --geogrid-file s3://jupiter-intern-projects/aarona/noah-mp-hue/milwaukee/geo_em.d01.milwaukee.nc
.PHONY: build clean push release run
