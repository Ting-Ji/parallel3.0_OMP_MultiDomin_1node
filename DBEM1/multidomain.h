#ifndef _MULTIDOMAIN_
#define _MULTIDOMAIN_

#include "dsquareelement.h"
#include "mat_vec.h"

#include <array>
#include <map>
#include <utility>
#include <vector>
#include <stdio.h>
#include <string>

struct Domain {
	int id;
	DynaMat mat;
	long eleBegin;
	long eleCount;
	long nodeBegin;
	long nodeCount;
	std::vector<long> elementIds;
	std::vector<long> boundaryNodeIds;
};

enum SurfaceType {
	SurfaceOuter = 0,
	SurfaceInterface = 1
};

struct InterfacePair {
	int domainA;
	int domainB;
	long eleA;
	long eleB;
	int normalSign;
	int localBForA[8];
	int localAForB[8];
	int nodeMapResolved;

	InterfacePair()
		: domainA(0), domainB(0), eleA(0), eleB(0), normalSign(-1), nodeMapResolved(0)
	{
		for (int i = 0; i < 8; ++i)
		{
			localBForA[i] = i;
			localAForB[i] = i;
		}
	}
};

struct DomainMaterialInput {
	int id;
	double E;
	double v;
	double Rou;
};

struct MultiDomainInputConfig {
	bool enabled;
	int domainCount;
	bool hasExplicitInterfaces;
	std::vector<DomainMaterialInput> materials;
	std::vector<InterfacePair> interfaces;

	MultiDomainInputConfig()
		: enabled(false), domainCount(1), hasExplicitInterfaces(false) {}
};

class DomainMaterialContext {
public:
	std::vector<DynaMat> mats;

	const DynaMat& Get(int domainId) const;
};

class MultiDomainModel {
public:
	bool enabled;
	std::vector<Domain> domains;
	std::vector<InterfacePair> interfaces;
	DomainMaterialContext materialContext;
	long elementCount;
	long nodeCount;

	MultiDomainModel();

	bool IsActive() const;
	int DomainCount() const;
	void PrintSummary(FILE* fp) const;
	bool Validate() const;
};

struct VariableRef {
	long globalBlock;
	double sign;
	bool known;

	VariableRef() : globalBlock(-1), sign(1.0), known(true) {}
	VariableRef(long block, double s, bool isKnown)
		: globalBlock(block), sign(s), known(isKnown) {}
};

class GlobalDofMap {
public:
	long unknownBlockCount;
	long knownBlockCount;
	long equationBlockCount;
	long outerDisplacementUnknownBlocks;
	long outerTractionUnknownBlocks;
	long interfaceDisplacementUnknownBlocks;
	long interfaceTractionUnknownBlocks;
	long historyUBlockCount;
	long historyTBlockCount;

	GlobalDofMap();

	bool Build(const MultiDomainModel& model, const DSquareElement* elements, long eleNum);
	VariableRef GetU(int domainId, long localNodeOrEle) const;
	VariableRef GetT(int domainId, long localNodeOrEle) const;
	VariableRef GetHistoryU(int domainId, long localNodeOrEle) const;
	VariableRef GetHistoryT(int domainId, long localNodeOrEle) const;
	long GetEquationRow(int domainId, long localNodeOrEle) const;
	long GetPreferredUnknownBlockForRow(long row) const;
	const char* GetUnknownDebugLabel(long unknownBlock) const;
	void PrintSummary(FILE* fp) const;
	bool Validate() const;

private:
	std::vector<VariableRef> m_uRefs;
	std::vector<VariableRef> m_tRefs;
	std::vector<VariableRef> m_historyURefs;
	std::vector<VariableRef> m_historyTRefs;
	std::vector<long> m_equationRows;
	std::vector<long> m_preferredUnknownBlockForRow;
	std::vector<std::string> m_unknownDebugLabels;
	std::vector<int> m_domainByElement;
};

class MultiDomainCCSRBuilder {
public:
	explicit MultiDomainCCSRBuilder(long blockRows = 0, long blockCols = 0);

	void Reset(long blockRows, long blockCols);
	void AddBlock(long rowBlock, long colBlock, double sign, const double block9[9]);
	void AddRhsBlock(long rowBlock, double sign, const double value3[3]);
	void MergeFrom(const MultiDomainCCSRBuilder& other);
	void Build(CCSRMat& out) const;
	void BuildSym(SymCCSRMat& out, double tolerance, bool* symmetric) const;
	void BuildRhs(Wvector& out) const;
	void BuildDense(std::vector<double>& out) const;
	long NonZeroBlocks() const;
	long BlockRows() const { return m_blockRows; }
	long BlockCols() const { return m_blockCols; }

private:
	long m_blockRows;
	long m_blockCols;
	std::map<std::pair<long, long>, std::array<double, 9> > m_blocks;
	std::vector<std::array<double, 3> > m_rhs;
};

struct MultiDomainState {
	Wvector globalUnknown;
	Wvector globalKnown;
	Wvector globalU;
	Wvector globalT;
	std::vector<Wvector> domainU;
	std::vector<Wvector> domainT;
};

int ReadOptionalMultiDomainInput(FILE* input, MultiDomainInputConfig& config);
int BuildMultiDomainModel(const MultiDomainInputConfig& config,
	long eleNum,
	long nodeNum,
	double defaultE,
	double defaultV,
	double defaultRou,
	double dt,
	MultiDomainModel& model);
void ApplyMultiDomainModelToElements(const MultiDomainModel& model, DSquareElement* elements, long eleNum);
int ResolveMultiDomainInterfaceNodeMaps(MultiDomainModel& model, const Point* nodeList, long nodeNum);

int AssembleMultiDomainStep0(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid,
	MultiDomainCCSRBuilder& unknownBuilder,
	MultiDomainCCSRBuilder& knownBuilder);

int AssembleMultiDomainHistoryStep(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long step,
	MultiDomainCCSRBuilder& mtsBuilder,
	MultiDomainCCSRBuilder& mgsBuilder);

int RunMultiDomainMatrixSelfCheck(DSquareElement* elements,
	const MultiDomainModel& model,
	const GlobalDofMap& dofMap,
	long** infElePid,
	long** elePid);

void WriteMultiDomainValidationMetrics(FILE* fp);

int DynaGaussSolverMultiDomainDense(DSquareElement* elements,
	BoundaryValue* bd,
	const MultiDomainModel& model,
	long NStep,
	double MaxLength,
	long** infElePid,
	long** elePid);

int DynaGMRESSolverMultiDomainCCSR(DSquareElement* elements,
	BoundaryValue* bd,
	const MultiDomainModel& model,
	long iterations,
	double error,
	long maxleafpointnum,
	long NStep,
	double MaxLength,
	long** infElePid,
	long** elePid);

#endif
