TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment

APP_SERVICE = app_service

APP_CXX = kmers-call-functions kmers-build-signatures
BIN_CXX = $(addprefix $(BIN_DIR)/,$(APP_CXX))
DEPLOY_CXX = $(addprefix $(TARGET)/bin,$(APP_CXX))


SRC_PERL = $(wildcard scripts/*.pl)
BIN_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_PERL))))
DEPLOY_PERL = $(addprefix $(TARGET)/bin/,$(basename $(notdir $(SRC_PERL))))

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)
PROD_TESTS = $(wildcard t/prod-tests/*.t)

STARMAN_WORKERS = 8
STARMAN_MAX_REQUESTS = 100

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) --define kb_service_dir=$(SERVICE_DIR) \
	--define kb_sphinx_port=$(SPHINX_PORT) --define kb_sphinx_host=$(SPHINX_HOST) \
	--define kb_starman_workers=$(STARMAN_WORKERS) \
	--define kb_starman_max_requests=$(STARMAN_MAX_REQUESTS)

all: bin

NuDB:
	git clone https://github.com/CPPAlliance/NuDB.git

bin: NuDB $(BIN_PERL) $(BIN_SERVICE_PERL) $(BIN_CXX)

#PROFILE = -pg
OPT = -O3
DEBUG = -g
INC = $(BOOST_INC) $(TBB_FLAGS) $(NUDB_INCLUDE) $(CMPH_INCLUDE)


CXXFLAGS = $(PROFILE) $(DEBUG) $(OPT) $(INC)
LDFLAGS = -Wl,-rpath,$(BOOST)/lib -Wl,-rpath,$(CMPH)/lib $(PROFILE)

LIBS = $(BOOST_LIBS) $(TBB_LIBS) $(CMPH_LIB)

BOOST = $(KB_RUNTIME)/boost-latest

BOOST_INC = -I$(BOOST)/include
BOOST_LIBS = \
	-L $(BOOST)/lib \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_regex

#
# Parallel algorithms on ubuntu result in TBB deprecation messages
#
TBB_FLAGS = -DTBB_SUPPRESS_DEPRECATED_MESSAGES=1 -DTBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS=1
TBB_LIBS = -ltbbmalloc -ltbb

CMPH = $(shell pwd)
CMPH_INCLUDE = -I$(CMPH)/include
CMPH_LIB = -L$(CMPH)/lib -lcmph

NUDB = NuDB
NUDB_INCLUDE = -I$(NUDB)/include

KMERS_CALL_FUNCTIONS_OBJS = src/kmers-call-functions.o src/fasta_parser.o
kmers-call-functions: NuDB $(KMERS_CALL_FUNCTIONS_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(KMERS_CALL_FUNCTIONS_OBJS) $(LIBS)

KMERS_MATRIX_DISTANCE_OBJS = src/kmers-matrix-distance.o src/fasta_parser.o
kmers-matrix-distance: $(KMERS_MATRIX_DISTANCE_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(KMERS_MATRIX_DISTANCE_OBJS) $(LIBS)

KMERS_BUILD_SIGNATURES = src/kmers-build-signatures.o src/fasta_parser.o
kmers-build-signatures: NuDB $(KMERS_BUILD_SIGNATURES)
	$(CXX) $(LDFLAGS) -o $@ $(KMERS_BUILD_SIGNATURES) $(LIBS)

tst-cmph: src/tst-cmph.o
	$(CXX) $(LDFLAGS) -o $@ src/tst-cmph.o $(LIBS)

write-cmph-from-kmers: src/write-cmph-from-kmers.o
	$(CXX) $(LDFLAGS) -o $@ src/write-cmph-from-kmers.o $(LIBS)

deploy: deploy-all
deploy-all: deploy-client 
deploy-client: deploy-libs deploy-scripts deploy-docs

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs: 


clean:


$(BIN_DIR)/%: service-scripts/%.pl $(TOP_DIR)/user-env.sh
	$(WRAP_PERL_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.py $(TOP_DIR)/user-env.sh
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

include $(TOP_DIR)/tools/Makefile.common.rules
