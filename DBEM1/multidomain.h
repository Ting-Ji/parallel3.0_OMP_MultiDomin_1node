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

// Internal zero-based contiguous range. Input V2 converts its 1-based start
// exactly once; legacy input resolves its ranges using the original split rule.
struct DomainElementRange {
	long begin;
	long count;
	DomainElementRange(long first = 0, long size = 0) : begin(first), count(size) {}
};

struct MultiDomainInputConfig {
	int formatVersion;
	bool enabled;
	int domainCount;
	bool hasExplicitInterfaces;
	std::vector<DomainMaterialInput> materials;
	std::vector<InterfacePair> interfaces;
	std::vector<DomainElementRange> elementRanges;

	MultiDomainInputConfig()
		: formatVersion(1), enabled(false), domainCount(1), hasExplicitInterfaces(false) {}
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

// Geometry-only snapshot. Rebuild after changing pairing, node IDs or local bases.
// A solve shares this immutable cache across assembly and state updates.
class MultiDomainInterfaceCache {
public:
	struct Entry {
		int domainId = -1;
		long interfaceId = -1;
		bool referenceSide = true;
		long node = -1, peerElement = -1, peerNode = -1;
		long referenceElement = -1, referenceNode = -1;
		std::array<double, 9> transform = {{1,0,0,0,1,0,0,0,1}};
	};
	MultiDomainInterfaceCache(const MultiDomainModel& model, const DSquareElement* elements);
	const Entry& Get(long element) const;
private:
	std::vector<Entry> entries;
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
	enum StorageMode { GlobalMap, IndependentRows };
	typedef std::map<long, std::array<double, 9> > BlockRow;
	// Matrix-only builders omit RHS storage until an RHS contribution is added.
	explicit MultiDomainCCSRBuilder(long blockRows = 0, long blockCols = 0, bool allocateRhs = true,
		StorageMode storage = GlobalMap);

	void Reset(long blockRows, long blockCols, bool allocateRhs = true, StorageMode storage = GlobalMap);
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
	bool HasIndependentRows() const { return m_storage == IndependentRows; }
	const BlockRow& Row(long row) const { return m_rows[(size_t)row]; }
	void ClearRows(long begin, long end);
	const std::array<double, 9>* FindBlock(long row, long col) const;
	template<class Visitor> void ForEachBlock(Visitor visit) const
	{
		if (HasIndependentRows())
		{
			for (long row = 0; row < m_blockRows; ++row)
				for (const auto& block : m_rows[(size_t)row])
					visit(row, block.first, block.second);
		}
		else
			for (const auto& block : m_blocks)
				visit(block.first.first, block.first.second, block.second);
	}

private:
	long m_blockRows;
	long m_blockCols;
	StorageMode m_storage;
	std::vector<BlockRow> m_rows;
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
