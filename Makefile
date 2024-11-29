include project.env
export

ifeq ($(IMG_NAME),)
	$(error IMG_NAME not set)
	exit 1
endif

$(info using image name: IMG_NAME=$(IMG_NAME))

REPO?=harbor.ngs-ai.com/docker
$(info using repo: REPO=$(REPO))

# decide the build command
BUILDER?=docker

ifeq ($(BUILDER),docker)
  $(info using docker (BUILDER=$(BUILDER)))
  BUILD_CMD := docker build
  PUSH_CMD := docker push
  TAG_CMD := docker tag
else ifeq ($(BUILDER),buildah)
  $(info using buildah (BUILDER=$(BUILDER)))
  # https://github.com/containers/buildah/issues/3790
  BUILD_CMD := buildah bud --isolation chroot
  PUSH_CMD := buildah push
  TAG_CMD := buildah tag
else
  $(error BUILDER=$(BUILDER) not supported)
  exit 1
endif

# get version
VERSION := $(shell git describe --exact-match HEAD)
VERSION_STATUS := $(shell git describe --tags --dirty --always)

# Build the container
build:
	$(BUILD_CMD) -t $(IMG_NAME) .

release: build publish 

publish: publish-latest publish-version

# Publish the `latest` taged container to ECR
publish-latest: tag-latest
	@echo 'publish latest to $(REPO)'
	$(PUSH_CMD) $(REPO)/$(IMG_NAME):latest

# Publish the `{version}` taged container to ECR
publish-version: tag-version 
	@echo 'publish $(VERSION) to $(REPO)'
	$(PUSH_CMD) $(REPO)/$(IMG_NAME):$(VERSION)

# Generate container tags for the `{version}` ans `latest` tags
tag: tag-latest tag-version

# Generate container `latest` tag
tag-latest:
	@echo 'create tag latest'
	$(TAG_CMD) $(IMG_NAME) $(REPO)/$(IMG_NAME):latest

# Generate container `version` tag, but check whether we're on a clean and tagged commit
tag-version:
	$(info version established as VERSION='$(VERSION)')
	$(info version status as VERSION_STATUS='$(VERSION_STATUS)')
	@if [ -z "$(VERSION)" ]; then\
        echo "No tagged version: you are not on a git tag. Exiting.";\
		exit 1;\
    fi
	@if echo "$(VERSION_STATUS)" | grep -q "dirty"; then\
        echo "You are on a dirty repository. Exiting.";\
		exit 1;\
    fi
	@echo 'create tag $(VERSION)'
	$(TAG_CMD) $(IMG_NAME) $(REPO)/$(IMG_NAME):$(VERSION)

# Below is used to check containers via trivy in the build process, via dumping to OCI
check-buildah:
	@if ! echo "$(BUILDER)" | grep -q "buildah"; then\
        echo "Not using BUILDER=buildah. Exiting.";\
		exit 1;\
    fi

dump: check-buildah
	buildah push $(IMG_NAME) oci-archive:./build.oci:$(IMG_NAME):latest

load: check-buildah
	buildah pull oci-archive:./build.oci:$(IMG_NAME):latest
