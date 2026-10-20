static void TestInterfaceCache()
{
    MultiDomainModel model;
    model.elementCount = model.nodeCount = 3;
    Domain da = {}, db = {};
    da.id = 0; da.elementIds = {0, 2};
    db.id = 1; db.elementIds = {1};
    model.domains = {da, db};
    InterfacePair pair;
    pair.domainA = 0; pair.domainB = 1;
    pair.eleA = 0; pair.eleB = 1; pair.nodeMapResolved = 1;
    model.interfaces.push_back(pair);
    DSquareElement elements[3];
    for (long i = 0; i < 3; ++i) elements[i].m_nodeID[0] = i;
    double a[9] = {0.6,0.8,0,-0.8,0.6,0,0,0,1};
    double b[9] = {1,0,0,0,0.8,0.6,0,-0.6,0.8};
    double* bases[3] = {a,b,nullptr};
    struct RestoreBases {
        double** previous;
        ~RestoreBases() { DSquareElement::m_transmat = previous; }
    } restore = {DSquareElement::m_transmat};
    DSquareElement::m_transmat = bases;
    const MultiDomainInterfaceCache cache(model, elements);
    CHECK(cache.Get(0).referenceSide && !cache.Get(1).referenceSide);
    CHECK(cache.Get(0).peerNode == 1 && cache.Get(1).peerNode == 0);
    CHECK(cache.Get(1).referenceElement == 0 && cache.Get(1).referenceNode == 0);
    CHECK(cache.Get(0).domainId == 0 && cache.Get(1).domainId == 1);
    CHECK(cache.Get(2).interfaceId == -1 && cache.Get(2).referenceSide);
    double expected[9];
    BuildLocalTransformFromReference(b, a, expected);
    for (int i = 0; i < 9; ++i) {
        CHECK(cache.Get(1).transform[(size_t)i] == expected[i]);
        CHECK(cache.Get(0).transform[(size_t)i] == (i % 4 == 0 ? 1.0 : 0.0));
        CHECK(cache.Get(2).transform[(size_t)i] == (i % 4 == 0 ? 1.0 : 0.0));
    }
    // Independent known product verifies orientation, not just the same helper.
    const double product[9] = {0.6,-0.8,0,0.64,0.48,0.6,-0.48,-0.36,0.8};
    for (int i = 0; i < 9; ++i) CHECK(fabs(expected[i] - product[i]) < 1e-14);
    DSquareElement::m_transmat = nullptr;
    CHECK(cache.Get(1).transform[0] == expected[0]); // Snapshot owns values.
    bool rejected = false;
    try { MultiDomainInterfaceCache invalid(model, elements); }
    catch (const std::runtime_error&) { rejected = true; }
    CHECK(rejected);
    DSquareElement::m_transmat = bases;
    model.interfaces.push_back(pair);
    rejected = false;
    try { MultiDomainInterfaceCache invalid(model, elements); }
    catch (const std::runtime_error&) { rejected = true; }
    CHECK(rejected);
    model.interfaces.clear();
    DSquareElement::m_transmat = nullptr;
    const MultiDomainInterfaceCache outerOnly(model, elements);
    CHECK(outerOnly.Get(1).interfaceId == -1);
}
