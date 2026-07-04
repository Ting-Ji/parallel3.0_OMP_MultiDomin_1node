#include "PthreadCCSR.h"
#ifdef _WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <string.h>

static long GmresDiagThreadCount = 0;
static int* GmresDiagSeen = 0;
static int* GmresDiagGroup = 0;
static int* GmresDiagProcessor = 0;

void GmresThreadDiagPrintSystemInfo(long requestedThreadNum)
{
	printf("GMRES requested thread_num = %ld\n", requestedThreadNum);
#ifdef _WIN32
	WORD groupCount = GetActiveProcessorGroupCount();
	printf("GMRES Windows processor groups = %u\n", (unsigned int)groupCount);
	for (WORD g = 0; g < groupCount; ++g)
		printf("GMRES group %u active processors = %lu\n", (unsigned int)g, (unsigned long)GetActiveProcessorCount(g));

	USHORT processGroupCount = 64;
	USHORT processGroups[64];
	memset(processGroups, 0, sizeof(processGroups));
	if (GetProcessGroupAffinity(GetCurrentProcess(), &processGroupCount, processGroups))
	{
		printf("GMRES process group affinity count = %u", (unsigned int)processGroupCount);
		for (USHORT i = 0; i < processGroupCount; ++i)
			printf(" %u", (unsigned int)processGroups[i]);
		printf("\n");
	}
	else
		printf("GMRES process group affinity unavailable, error = %lu\n", (unsigned long)GetLastError());

	DWORD_PTR processMask = 0;
	DWORD_PTR systemMask = 0;
	if (GetProcessAffinityMask(GetCurrentProcess(), &processMask, &systemMask))
		printf("GMRES process affinity mask = 0x%p, system affinity mask = 0x%p\n", (void*)processMask, (void*)systemMask);
	else
		printf("GMRES process affinity mask unavailable, error = %lu\n", (unsigned long)GetLastError());
#else
	printf("GMRES processor group diagnostics are only available on Windows.\n");
#endif
}

void GmresThreadDiagBegin(long requestedThreadNum)
{
	delete[] GmresDiagSeen;
	delete[] GmresDiagGroup;
	delete[] GmresDiagProcessor;
	GmresDiagSeen = 0;
	GmresDiagGroup = 0;
	GmresDiagProcessor = 0;
	GmresDiagThreadCount = requestedThreadNum;
	if (requestedThreadNum <= 0)
		return;
	GmresDiagSeen = new int[requestedThreadNum];
	GmresDiagGroup = new int[requestedThreadNum];
	GmresDiagProcessor = new int[requestedThreadNum];
	for (long i = 0; i < requestedThreadNum; ++i)
	{
		GmresDiagSeen[i] = 0;
		GmresDiagGroup[i] = -1;
		GmresDiagProcessor[i] = -1;
	}
}

void GmresThreadDiagRecord(long threadID)
{
	if (threadID < 0 || threadID >= GmresDiagThreadCount || !GmresDiagSeen)
		return;
#ifdef _WIN32
	PROCESSOR_NUMBER proc;
	GetCurrentProcessorNumberEx(&proc);
	GmresDiagGroup[threadID] = (int)proc.Group;
	GmresDiagProcessor[threadID] = (int)proc.Number;
#else
	GmresDiagGroup[threadID] = -1;
	GmresDiagProcessor[threadID] = -1;
#endif
	GmresDiagSeen[threadID] = 1;
}

void GmresThreadDiagPrintCreateResult(long requestedThreadNum, long successCount, const int* returnCodes)
{
	printf("GMRES pthread create requested = %ld, success = %ld\n", requestedThreadNum, successCount);
	if (!returnCodes)
		return;
	for (long i = 0; i < requestedThreadNum; ++i)
	{
		if (returnCodes[i] != 0)
			printf("GMRES pthread create failed: thread_id = %ld, return_code = %d\n", i, returnCodes[i]);
	}
}

void GmresThreadDiagEnd()
{
	if (!GmresDiagSeen || GmresDiagThreadCount <= 0)
		return;

	long recorded = 0;
	for (long i = 0; i < GmresDiagThreadCount; ++i)
	{
		if (GmresDiagSeen[i])
			++recorded;
	}
	printf("GMRES matrix worker recorded threads = %ld / %ld\n", recorded, GmresDiagThreadCount);

#ifdef _WIN32
	WORD activeGroups = GetActiveProcessorGroupCount();
	for (WORD g = 0; g < activeGroups; ++g)
	{
		long count = 0;
		for (long i = 0; i < GmresDiagThreadCount; ++i)
		{
			if (GmresDiagSeen[i] && GmresDiagGroup[i] == (int)g)
				++count;
		}
		printf("GMRES worker group %u threads = %ld\n", (unsigned int)g, count);
	}

	int uniqueCount = 0;
	int uniquePairs[4096];
	for (long i = 0; i < GmresDiagThreadCount; ++i)
	{
		if (!GmresDiagSeen[i])
			continue;
		int key = GmresDiagGroup[i] * 1024 + GmresDiagProcessor[i];
		int exists = 0;
		for (int j = 0; j < uniqueCount; ++j)
		{
			if (uniquePairs[j] == key)
			{
				exists = 1;
				break;
			}
		}
		if (!exists && uniqueCount < 4096)
			uniquePairs[uniqueCount++] = key;
	}
	printf("GMRES unique logical processors touched = %d\n", uniqueCount);

	for (long i = 0; i < GmresDiagThreadCount; ++i)
	{
		if (GmresDiagSeen[i])
			printf("GMRES worker thread %ld ran on group %d processor %d\n", i, GmresDiagGroup[i], GmresDiagProcessor[i]);
		else
			printf("GMRES worker thread %ld did not record processor location\n", i);
	}
#endif
}





void *PthGetMatrixDynaUnknown0ij(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	GmresThreadDiagRecord(ID);
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest= NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j, m_InfElePID, m_ElePID,ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MCCSRtemp.assign(AssistCoef[ID][l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 0;
}

void* PthGetMatrixDynaUnknown0ijT0(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}



	//printf("�̺߳�:%ld\n", ID);

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED

	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i ;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = UnknownSubMatrix(m_DSE, i, j, m_InfElePID, m_ElePID, ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MT[0].assign(AssistCoef[ID][l], pos+1, 3 * k+1, 3, 3);
				}
			}
		}
	}


	return 0;
}

void *PthGetMatrixDynaKnown0ij(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i ;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = knownSubMatrix(m_DSE,i, j, m_InfElePID, m_ElePID, ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MCCSRtemp.assign(AssistCoef[ID][l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	return 0;
}

void *PthGetMatrixDynaT1ij(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos, ii, jj;;
	int Flag;
	int ipos;
	double temp[8][9];
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			for (ii = 0; ii < 8; ++ii) {
				for (jj = 0; jj < 9; ++jj) {
					temp[ii][jj] = 0.0;
				}
			}

			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

					for (jj = 0; jj < 9; ++jj) {
						temp[l][jj] += AssistTij[ID][l][jj];
					}
				}
			}

			////���޵�Ԫ����
			//m_InfEleFlag = m_DSE[j].InfEleFlag;
			//if (m_InfEleFlag != (-1))
			//{
			//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
			//	}
			//	//����
			//	for (l = 0; l < 1; ++l)
			//	{
			//		k = m_DSE[j].m_nodeID[l];
			//		// sub = DSquareElement::m_AssistTij[l];
			//		for (jj = 0; jj < 9; ++jj) {
			//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];

			//		}
			//		//system("pause");
			//		//A.assign(temp[l], pos, 3 * k, 3, 3);
			//		//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

			//	}
			//}

			for (l = 0; l < 1; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
			}
		}

	}
	return 0;
}

void* PthGetMatrixDynaT1ijNew(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;

	long i, j, k, l, pos, ii, jj;;
	int Flag;
	double temp[8][9];

	for (i = start; i < end; ++i)
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			for (ii = 0; ii < 8; ++ii) {
				for (jj = 0; jj < 9; ++jj) {
					temp[ii][jj] = 0.0;
				}
			}

			Flag = m_DSE[j].IntDynaTijJudge(1, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistTij[l];
					//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

					for (jj = 0; jj < 9; ++jj) {
						temp[l][jj] += AssistTij[ID][l][jj];
					}
				}
			}

			////���޵�Ԫ����
			//m_InfEleFlag = m_DSE[j].InfEleFlag;
			//if (m_InfEleFlag != (-1))
			//{
			//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
			//	}
			//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
			//		m_DSE[j].IntDynaTijPWforInfEle(1, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
			//	}
			//	//����
			//	for (l = 0; l < 1; ++l)
			//	{
			//		k = m_DSE[j].m_nodeID[l];
			//		// sub = DSquareElement::m_AssistTij[l];
			//		for (jj = 0; jj < 9; ++jj) {
			//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];

			//		}
			//		//system("pause");
			//		//A.assign(temp[l], pos, 3 * k, 3, 3);
			//		//A.assign(DSquareElement::m_AssistTij[l], pos + 1, 3 * k + 1, 3, 3);

			//	}
			//}

			for (l = 0; l < 1; ++l)
			{
				k = m_DSE[j].m_nodeID[l];
				MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
			}
		}

	}
	return 0;
}

void *PthGetMatrixDynaT2ij(void* arg)
{

	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos, ii, jj;;
	int Flag;
	int ipos;
	double temp[8][9];
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}


	//printf("�̺߳�:%ld\n", ID);

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED

	if (StepNum > 1.01) // �������
	{
		for (i = count; i < end; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 8; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}

				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += AssistTij[ID][l][jj];
						}
					}
				}
				////���޵�Ԫ����
				//m_InfEleFlag = m_DSE[j].InfEleFlag;
				//if (m_InfEleFlag != (-1))
				//{
				//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				//	}
				//	//����
				//	for (l = 0; l < 1; ++l)
				//	{
				//		k = m_DSE[j].m_nodeID[l];
				//		// sub = DSquareElement::m_AssistTij[l];
				//		for (jj = 0; jj < 9; ++jj) {
				//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
				//		}
				//		//A.assign(temp[l], pos, 3 * k, 3, 3);
				//	}
				//	//system("pause");
				//}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	else // n = 1
	{
		for (i = count; i < end; ++i)//recycle for source element
		{
			pos = 3 * i ;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 8; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 1; ++l)
					{
						AssistTij[ID][l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						AssistTij[ID][l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						AssistTij[ID][l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						AssistTij[ID][l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						AssistTij[ID][l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						AssistTij[ID][l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						AssistTij[ID][l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						AssistTij[ID][l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						AssistTij[ID][l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += AssistTij[ID][l][jj];
						}
					}
				}
				////���޵�Ԫ����
				//m_InfEleFlag = m_DSE[j].InfEleFlag;
				//if (m_InfEleFlag != (-1))
				//{
				//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				//	}
				//	//����
				//	for (l = 0; l < 1; ++l)
				//	{
				//		k = m_DSE[j].m_nodeID[l];
				//		// sub = DSquareElement::m_AssistTij[l];
				//		for (jj = 0; jj < 9; ++jj) {
				//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
				//		}
				//		//A.assign(DSquareElement::m_AssistTij[l], pos , 3 * k , 3, 3);
				//	}
				//}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
				}


			}
		}
	}
	return 0;
}

void* PthGetMatrixDynaT2ijNew(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;
	int ipos;
	long i, j, k, l, pos, ii, jj;;
	int Flag;
	double temp[8][9];

	if (StepNum > 1.01) // �������
	{
		for (i = start; i < end; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 8; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}

				Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += AssistTij[ID][l][jj];
						}
					}
				}
				////���޵�Ԫ����
				//m_InfEleFlag = m_DSE[j].InfEleFlag;
				//if (m_InfEleFlag != (-1))
				//{
				//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				//	}
				//	//����
				//	for (l = 0; l < 1; ++l)
				//	{
				//		k = m_DSE[j].m_nodeID[l];
				//		// sub = DSquareElement::m_AssistTij[l];
				//		for (jj = 0; jj < 9; ++jj) {
				//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
				//		}
				//		//A.assign(temp[l], pos, 3 * k, 3, 3);
				//	}
				//	//system("pause");
				//}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	else // n = 1
	{
		for (i = start; i < end; ++i)//recycle for source element
		{
			pos = 3 * i;
			for (j = 0; j < EleNum; ++j)
			{
				for (ii = 0; ii < 8; ++ii) {
					for (jj = 0; jj < 9; ++jj) {
						temp[ii][jj] = 0.0;
					}
				}
				if (m_DSE[j].IsIn(i, ipos))
				{
					for (l = 0; l < 1; ++l)
					{
						AssistTij[ID][l][0] = m_DSE[j].m_OnElementT2ij[ipos][l][0];
						AssistTij[ID][l][1] = m_DSE[j].m_OnElementT2ij[ipos][l][1];
						AssistTij[ID][l][2] = m_DSE[j].m_OnElementT2ij[ipos][l][2];
						AssistTij[ID][l][3] = m_DSE[j].m_OnElementT2ij[ipos][l][3];
						AssistTij[ID][l][4] = m_DSE[j].m_OnElementT2ij[ipos][l][4];
						AssistTij[ID][l][5] = m_DSE[j].m_OnElementT2ij[ipos][l][5];
						AssistTij[ID][l][6] = m_DSE[j].m_OnElementT2ij[ipos][l][6];
						AssistTij[ID][l][7] = m_DSE[j].m_OnElementT2ij[ipos][l][7];
						AssistTij[ID][l][8] = m_DSE[j].m_OnElementT2ij[ipos][l][8];
					}

					Flag = 1;
				}
				else
					Flag = m_DSE[j].IntDynaTijJudge(2, m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

				if (Flag)
				{
					for (l = 0; l < 1; ++l)
					{
						k = m_DSE[j].m_nodeID[l];
						for (jj = 0; jj < 9; ++jj) {
							temp[l][jj] += AssistTij[ID][l][jj];
						}
					}
				}
				////���޵�Ԫ����
				//m_InfEleFlag = m_DSE[j].InfEleFlag;
				//if (m_InfEleFlag != (-1))
				//{
				//	if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][5]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 1);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][6]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 2);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][7]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 3);
				//	}
				//	else if (InfElePID[m_InfEleFlag][2] == m_ElePID[j][8]) {
				//		m_DSE[j].IntDynaTijPWforInfEle(2, m_DSE[j].m_nodelist[i], n, DSquareElement::m_DMat.Dt, 4);
				//	}
				//	//����
				//	for (l = 0; l < 1; ++l)
				//	{
				//		k = m_DSE[j].m_nodeID[l];
				//		// sub = DSquareElement::m_AssistTij[l];
				//		for (jj = 0; jj < 9; ++jj) {
				//			temp[l][jj] += DSquareElement::m_AssistTij[l][jj];
				//		}
				//		//A.assign(DSquareElement::m_AssistTij[l], pos , 3 * k , 3, 3);
				//	}
				//}

				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					MCCSRtemp.assign(temp[l], pos, 3 * k, 3, 3);
				}


			}
		}
	}
	return 0;
}

void *PthGetMatrixDynaUij(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}


	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED

	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], StepNum, dt_p,ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MCCSRtemp.assign(AssistUij[ID][l], pos, 3 * k, 3, 3);
				}
			}
		}
	}

	return 0;
}

void* PthGetMatrixDynaSymUij(void* arg)
{

	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;


	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED

	
	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);


			/*if (Flag == 0) {
				FlagUij[StepNum][i][0]++;
			}
			else if (Flag == 1) {
				FlagUij[StepNum][i][1]++;
			}
			else if (Flag == 2) {
				FlagUij[StepNum][i][2]++;
			}*/


			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MCCSRtemp.assignSym(AssistUij[ID][l], pos, 3 * k, 3, 3);
					//MSymCtemp.assign(AssistUij[ID][l], pos, 3 * k, 3, 3);
				}
			}
		}



	}




	return 0;
}

void* PthGetMatrixDynaSymUijNew(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;

	long i, j, k, l, pos, ii, jj;;
	int Flag;
	double temp[8][9];

	for (i = start; i < end; ++i)//recycle for source element   i->source ID
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaUijJudge(m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);

			if (Flag)
			{
				for (l = 0; l < 1; ++l)
				{
					k = m_DSE[j].m_nodeID[l];
					// sub = DSquareElement::m_AssistCoef[l];
					MCCSRtemp.assignSym(AssistUij[ID][l], pos, 3 * k, 3, 3);
					//MSymCtemp.assign(AssistUij[ID][l], pos, 3 * k, 3, 3);
				}
			}
		}
	}
	return 0;
}


void* PthSparseConvertCCSR(void* arg)
{
	PthArg* PthTemp= (PthArg*)arg;
	long ID = PthTemp->targs;
	long start= PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;
	
	for (i = start; i < end; ++i)
	{
			for (j = 0; j < MCCSRtemp.row[i] * 9; ++j) {
				MTS[StepNum].value[nnum] = MCCSRtemp.a[i][j];
				++nnum;
			}
			for (j = 0; j < MCCSRtemp.row[i]; j++) {
				MTS[StepNum].col_index[Pcolnum] = MCCSRtemp.col[i][j];
				Pcolnum++;
			}
			sum += MCCSRtemp.row[i];
			MTS[StepNum].row_ptr[i + 1] = sum ;
	}
	for (i = start; i < end; ++i)
	{
			for (int j = 0; j < MCCSRtemp.n; ++j) {
				MCCSRtemp.col[i][j] = 0;
			}
			for (int j = 0; j < MCCSRtemp.n * 9; ++j) {
				MCCSRtemp.a[i][j] = 0;
			}
			MCCSRtemp.row[i] = 0;
	}
	if (ID = thread_num - 1)
		MCCSRtemp.num = 0;

	return 0;
}

void* PthSparseAddConvertCCSR(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->AddStrat;
	long end = PthTemp->AddEnd;
	long i, j;
	long nnum = (PthTemp->AddSum) * 9;
	long Pcolnum = PthTemp->AddSum;
	long sum = Pcolnum;

	for (i = start; i < end; ++i)
	{
		for (j = 0; j < MCCSRtemp.row[i] * 9; ++j) {
			MTS[StepNum].value[nnum] = MCCSRtemp.a[i][j];
			++nnum;
		}
		for (j = 0; j < MCCSRtemp.row[i]; j++) {
			MTS[StepNum].col_index[Pcolnum] = MCCSRtemp.col[i][j];
			Pcolnum++;
		}
		sum += MCCSRtemp.row[i];
		MTS[StepNum].row_ptr[i + 1] = sum;
	}
	for (i = start; i < end; ++i)
	{
		for (int j = 0; j < MCCSRtemp.n; ++j) {
			MCCSRtemp.col[i][j] = 0;
		}
		for (int j = 0; j < MCCSRtemp.n * 9; ++j) {
			MCCSRtemp.a[i][j] = 0;
		}
		MCCSRtemp.row[i] = 0;
	}
	if (ID = thread_num - 1)
		MCCSRtemp.num = 0;

	return 0;
}


void* PthSparseConvertCCSRMTS0(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;

	for (i = start; i < end; ++i)
	{
		for (j = 0; j < MCCSRtemp.row[i] * 9; ++j) {
			MTS0.value[nnum] = MCCSRtemp.a[i][j];
			++nnum;
		}
		for (j = 0; j < MCCSRtemp.row[i]; j++) {
			MTS0.col_index[Pcolnum] = MCCSRtemp.col[i][j];
			Pcolnum++;
		}
		sum += MCCSRtemp.row[i];
		MTS0.row_ptr[i + 1] = sum;
	}
	for (i = start; i < end; ++i)
	{
		for (int j = 0; j < MCCSRtemp.n; ++j) {
			MCCSRtemp.col[i][j] = 0;
		}
		for (int j = 0; j < MCCSRtemp.n * 9; ++j) {
			MCCSRtemp.a[i][j] = 0;
		}
		MCCSRtemp.row[i] = 0;
	}
	if (ID = thread_num)
		MCCSRtemp.num = 0;

	return 0;
}

void* PthSparseConvertSymCCSR(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	long nnum = (PthTemp->sum) * 6;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;

	for (i = start; i < end; ++i)
	{
		for (j = 0; j < MCCSRtemp.row[i] * 6; ++j) {
			MGS[StepNum].value[nnum] = MCCSRtemp.a[i][j];
			++nnum;
		}
		for (j = 0; j < MCCSRtemp.row[i]; j++) {
			MGS[StepNum].col_index[Pcolnum] = MCCSRtemp.col[i][j];
			Pcolnum++;
		}
		sum += MCCSRtemp.row[i];
		MGS[StepNum].row_ptr[i + 1] = sum;
	}
	for (i = start; i < end; ++i)
	{
		for (int j = 0; j < MCCSRtemp.n; ++j) {
			MCCSRtemp.col[i][j] = 0;
		}
		for (int j = 0; j < MCCSRtemp.n * 6; ++j) {
			MCCSRtemp.a[i][j] = 0;
		}
		MCCSRtemp.row[i] = 0;
	}
	if (ID = thread_num)
		MCCSRtemp.num = 0;

	return 0;
}


void* PthSparseAdditionCCSR(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j, k;
	long nnum = (PthTemp->sum) * 9;
	long Pcolnum = PthTemp->sum;
	long sum = Pcolnum;
	long b_row;//b�з���Ԫ�ظ���
	int A_index, M_index;//����M=A+Mһ���е�A��M��ָʾ��
	int M_row_ptr;
	int M_rownum;//M�е�i�з���Ԫ�ظ���

	WvectorCCSR b;
	b.initial(MCCSRtemp.m * 3);

	for (i = start; i < end; ++i)
	{
		A_index = 0;
		M_index = 0;
		M_row_ptr = MTS[StepNum].row_ptr[i];//��i�е�һ������Ԫ����M.value�е�λ�ã�
		b_row = 0;
		M_rownum = MTS[StepNum].row_ptr[i + 1] - MTS[StepNum].row_ptr[i];
		b.reinitial();//b��0
		while (A_index < MCCSRtemp.row[i] && M_index < M_rownum) {
			if (MCCSRtemp.col[i][A_index] < MTS[StepNum].col_index[M_row_ptr + M_index]) {//A<B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = MCCSRtemp.a[i][A_index * 9 + k];
				}

				b.col[b_row] = MCCSRtemp.col[i][A_index];
				++b_row;
				++A_index;

			}
			else if (MCCSRtemp.col[i][A_index] > MTS[StepNum].col_index[M_row_ptr + M_index]) {//A>B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = MTS[StepNum].value[(M_row_ptr + M_index) * 9 + k];
				}
				b.col[b_row] = MTS[StepNum].col_index[M_row_ptr + M_index];
				++b_row;
				++M_index;
			}
			else {//A=B
				for (k = 0; k < 9; k++) {
					b.value[b_row * 9 + k] = MCCSRtemp.a[i][A_index * 9 + k] + MTS[StepNum].value[(M_row_ptr + M_index) * 9 + k];
				}
				b.col[b_row] = MCCSRtemp.col[i][A_index];
				++b_row;
				++M_index;
				++A_index;
			}
		}
		while (A_index < MCCSRtemp.row[i]) {

			for (k = 0; k < 9; k++) {
				b.value[b_row * 9 + k] = MCCSRtemp.a[i][A_index * 9 + k];
			}
			b.col[b_row] = MCCSRtemp.col[i][A_index];
			++b_row;
			++A_index;
		}
		while (M_index < M_rownum) {
			for (k = 0; k < 9; k++) {
				b.value[b_row * 9 + k] = MTS[StepNum].value[(M_row_ptr + M_index) * 9 + k];
			}
			b.col[b_row] = MTS[StepNum].col_index[M_row_ptr + M_index];
			++b_row;
			++M_index;
		}
		MCCSRtemp.row[i] = b_row;
		for (j = 0; j < b_row; ++j) {
			for (k = 0; k < 9; k++) {
				MCCSRtemp.a[i][j * 9 + k] = b.value[j * 9 + k];
			}
			MCCSRtemp.col[i][j] = b.col[j];
		}
	}


	return 0;
}

void* PthDistr(int* row, PthArg* Thread,long num, long n)
{
	long averg;
	long sum=0;
	long flag = 0;
	long count = 0;
	for (long i = 0; i < thread_num - 1; i++){
		averg = (num - count) / (thread_num - i);
		Thread[i].start = flag;
		Thread[i].sum = count;
		sum = 0;
		while (sum < averg) {
			sum += row[flag];
			flag++;
		}
		count += sum;
		Thread[i].end = flag;
		sum = 0;
	}

	Thread[thread_num - 1].start = flag;
	Thread[thread_num - 1].end = n;
	Thread[thread_num - 1].sum = count;
	return 0;
}

void* PthDistrAdd(int* row, PthArg* Thread, long num, long n)
{
	/*long averg;
	long sum = 0;
	long flag = 0;
	long count = 0;
	for (long i = 0; i < thread_num - 1; i++) {
		averg = (num - count) / (thread_num - i);
		ThreadAdd[i].strat = flag;
		ThreadAdd[i].sum = count;
		sum = 0;
		while (sum < averg) {
			sum += row[flag];
			flag++;
		}
		count += sum;
		ThreadAdd[i].end = flag;
		sum = 0;
	}

	ThreadAdd[thread_num - 1].strat = flag;
	ThreadAdd[thread_num - 1].end = n;
	ThreadAdd[thread_num - 1].sum = count;
	return nullptr;*/



	long averg;
	long sum = 0;
	long flag = 0;
	long count = 0;
	for (long i = 0; i < thread_num - 1; i++) {
		averg = (num - count) / (thread_num - i);
		Thread[i].AddStrat = flag;
		Thread[i].AddSum = count;
		sum = 0;
		while (sum < averg) {
			sum += row[flag];
			flag++;
		}
		count += sum;
		Thread[i].AddEnd = flag;
		sum = 0;
	}

	Thread[thread_num - 1].AddStrat = flag;
	Thread[thread_num - 1].AddEnd = n;
	Thread[thread_num - 1].AddSum = count;
	return 0;
}

void* PthDistr(PthArg* Thread)
{
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;

	for (long i = 0; i < thread_num ; i++) {
		Thread[i].start = i * size;
		if (i < rest) {
			Thread[i].start += i;
			Thread[i].end = Thread[i].start + size + 1;
		}
		else {
			Thread[i].start += rest;
			Thread[i].end = Thread[i].start + size;
		}

	}
	return 0;

}



void* PthAllTraverse(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	//printf("�̺߳�:%ld\n", ID);
	long temp;

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	// TEST CODE: NECESSARY, OTHERWISE SOME ENTRIES OF A NOT ASSIGNED
	for (i = count; i < end; ++i)//recycle for source element   i->source ID
	{
		Pthnum[i] = 0;
		temp = 0;
		for (j = 0; j < EleNum; ++j)
		{
			Flag = m_DSE[j].IntDynaJudge(m_DSE[j].m_nodelist[i], StepNum, dt_p, ID);
			temp += Flag;
		}
		Pthnum[i] = temp;
	}

	return 0;
}

void* PthSpMVKnown(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	//int i, j, k;
	//int n = NodeNum;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (i = start; i < end; ++i) {
		Bvector[i * 3 + 1] = 0.0;
		Bvector[i * 3 + 2] = 0.0;
		Bvector[i * 3 + 3] = 0.0;
		rownum = MTS[0].row_ptr[i + 1] - MTS[0].row_ptr[i];
		rowstart = MTS[0].row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = bd[StepNum].lu[MTS[0].col_index[rowstart + j] * 3 + 1];
			btemp[1] = bd[StepNum].lu[MTS[0].col_index[rowstart + j] * 3 + 2];
			btemp[2] = bd[StepNum].lu[MTS[0].col_index[rowstart + j] * 3 + 3];
			Bvector[i * 3 + 1] += MTS[0].value[rowstart * 9 + j * 9] * btemp[0] + MTS[0].value[rowstart * 9 + j * 9 + 1] * btemp[1] + MTS[0].value[rowstart * 9 + j * 9 + 2] * btemp[2];
			Bvector[i * 3 + 2] += MTS[0].value[rowstart * 9 + j * 9 + 3] * btemp[0] + MTS[0].value[rowstart * 9 + j * 9 + 4] * btemp[1] + MTS[0].value[rowstart * 9 + j * 9 + 5] * btemp[2];
			Bvector[i * 3 + 3] += MTS[0].value[rowstart * 9 + j * 9 + 6] * btemp[0] + MTS[0].value[rowstart * 9 + j * 9 + 7] * btemp[1] + MTS[0].value[rowstart * 9 + j * 9 + 8] * btemp[2];
		}
	}

	return 0;
}

void* PthSpMVUnkown(void* arg)
{
	return 0;
}

void* PthSpMVT(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	//int i, j, k;
	//int n = NodeNum;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (i = start; i < end; ++i) {
		CB2[i * 3 + 1] = 0.0;
		CB2[i * 3 + 2] = 0.0;
		CB2[i * 3 + 3] = 0.0;
		rownum = MTS[ConvoluNum].row_ptr[i + 1] - MTS[ConvoluNum].row_ptr[i];//某细胞行非零细胞矩阵个数
		rowstart = MTS[ConvoluNum].row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = CB1[MTS[ConvoluNum].col_index[rowstart + j] * 3 + 1];
			btemp[1] = CB1[MTS[ConvoluNum].col_index[rowstart + j] * 3 + 2];
			btemp[2] = CB1[MTS[ConvoluNum].col_index[rowstart + j] * 3 + 3];
			CB2[i * 3 + 1] += MTS[ConvoluNum].value[rowstart * 9 + j * 9] * btemp[0] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 1] * btemp[1] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 2] * btemp[2];
			CB2[i * 3 + 2] += MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 3] * btemp[0] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 4] * btemp[1] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 5] * btemp[2];
			CB2[i * 3 + 3] += MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 6] * btemp[0] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 7] * btemp[1] + MTS[ConvoluNum].value[rowstart * 9 + j * 9 + 8] * btemp[2];
		}
	}

	return 0;
}

void* PthSpMVG(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, j;
	//int i, j, k;
	//int n = NodeNum;
	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (i = start; i < end; ++i) {
		CB2[i * 3 + 1] = 0.0;
		CB2[i * 3 + 2] = 0.0;
		CB2[i * 3 + 3] = 0.0;
		rownum = MGS[ConvoluNum].row_ptr[i + 1] - MGS[ConvoluNum].row_ptr[i];
		rowstart = MGS[ConvoluNum].row_ptr[i];
		for (j = 0; j < rownum; ++j) {
			btemp[0] = CB1[MGS[ConvoluNum].col_index[rowstart + j] * 3 + 1];
			btemp[1] = CB1[MGS[ConvoluNum].col_index[rowstart + j] * 3 + 2];
			btemp[2] = CB1[MGS[ConvoluNum].col_index[rowstart + j] * 3 + 3];
			CB2[i * 3 + 1] += MGS[ConvoluNum].value[rowstart * 6 + j * 6] * btemp[0] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 1] * btemp[1] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 2] * btemp[2];
			CB2[i * 3 + 2] += MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 1] * btemp[0] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 3] * btemp[1] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 4] * btemp[2];
			CB2[i * 3 + 3] += MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 2] * btemp[0] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 4] * btemp[1] + MGS[ConvoluNum].value[rowstart * 6 + j * 6 + 5] * btemp[2];
		}
	}

	return 0;

}

void* PthNorm(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = (*PthTemp).targs;
	long i, j, k, l, pos;
	int Flag;
	long size = NodeNum / thread_num;
	long rest = NodeNum % thread_num;
	long count = ID * size;
	long end;
	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}


	for (i = count; i < end; ++i)
	{
		pos = 3 * i;
		for (j = 0; j < EleNum; ++j)
		{

		}
	}

	return 0;
}

void* PthbdConvert(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, k;

	for (i = start; i < end; ++i) {
		k = i;//k = 8 * i;
		bd[StepNum].Convert(m_DSE[i].BCID, k + 0);
		//Convert(m_DSE[i].BCID, k + 1);
		//Convert(m_DSE[i].BCID, k + 2);
		//Convert(m_DSE[i].BCID, k + 3);
		//Convert(m_DSE[i].BCID, k + 4);
		//Convert(m_DSE[i].BCID, k + 5);
		//Convert(m_DSE[i].BCID, k + 6);
		//Convert(m_DSE[i].BCID, k + 7);
		
	}

	return 0;
}

void* PthbdAllLocToAbs(void* arg)
{
	PthArg* PthTemp = (PthArg*)arg;
	long ID = PthTemp->targs;
	long start = PthTemp->start;
	long end = PthTemp->end;
	long i, k;

	for (i = start; i < end; ++i) {
		bd[StepNum].LocToAbs(i, DSquareElement::m_transmat[i]);

	}

	return 0;
}



void* SpmvConvertConvol(PthArg** Spmvthread, PthConvol* Convolthread,long MaxN)
{
	long i, j;

	for (i = 0; i < thread_num; i++) {
		Convolthread[i].targs = i;
	}

	for (i = 0; i < MaxN + 1; ++i) {
		for (j = 0; j < thread_num; ++j) {
			Convolthread[j].start[i] = Spmvthread[i][j].start;
			Convolthread[j].end[i] = Spmvthread[i][j].end;
		}
	}

	return nullptr;
}

void* PthConvolutinT(void* arg)
{
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = PthTemp->targs;
	long* start = PthTemp->start;
	long* end = PthTemp->end;
	long MaxN = PthTemp->MaxN;
	long Cnum, i, j;
	//Wvector CBthread;
	//CBthread.initial(3 * NodeNum);
	//CBthread = 0;

	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (Cnum = 1; Cnum <= StepNum; ++Cnum) {
		if (Cnum <= MaxN) {
			for (i = start[Cnum]; i < end[Cnum]; ++i) {
				rownum = MTS[Cnum].row_ptr[i + 1] - MTS[Cnum].row_ptr[i];//某细胞行非零细胞矩阵个数
				rowstart = MTS[Cnum].row_ptr[i];
				for (j = 0; j < rownum; ++j) {
					btemp[0] = CBT[Cnum][MTS[Cnum].col_index[rowstart + j] * 3 + 1];
					btemp[1] = CBT[Cnum][MTS[Cnum].col_index[rowstart + j] * 3 + 2];
					btemp[2] = CBT[Cnum][MTS[Cnum].col_index[rowstart + j] * 3 + 3];
					Bvector[ID][i * 3 + 1] += MTS[Cnum].value[rowstart * 9 + j * 9] * btemp[0] + MTS[Cnum].value[rowstart * 9 + j * 9 + 1] * btemp[1] + MTS[Cnum].value[rowstart * 9 + j * 9 + 2] * btemp[2];
					Bvector[ID][i * 3 + 2] += MTS[Cnum].value[rowstart * 9 + j * 9 + 3] * btemp[0] + MTS[Cnum].value[rowstart * 9 + j * 9 + 4] * btemp[1] + MTS[Cnum].value[rowstart * 9 + j * 9 + 5] * btemp[2];
					Bvector[ID][i * 3 + 3] += MTS[Cnum].value[rowstart * 9 + j * 9 + 6] * btemp[0] + MTS[Cnum].value[rowstart * 9 + j * 9 + 7] * btemp[1] + MTS[Cnum].value[rowstart * 9 + j * 9 + 8] * btemp[2];

				}

			}
		
		}
		

	}
		return nullptr;
}



void* PthConvolutinG(void* arg)
{
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = PthTemp->targs;
	long* start = PthTemp->start;
	long* end = PthTemp->end;
	long MaxN = PthTemp->MaxN;
	long Cnum, i, j;
	//Wvector CBthread;
	//CBthread.initial(3 * NodeNum);
	//CBthread = 0;

	int rownum;
	int rowstart;
	double btemp[3];
	btemp[0] = btemp[1] = btemp[2] = 0;

	for (Cnum = 1; Cnum <= StepNum; ++Cnum) {
		if (Cnum <= MaxN) {
			for (i = start[Cnum]; i < end[Cnum]; ++i) {
				rownum = MGS[Cnum].row_ptr[i + 1] - MGS[Cnum].row_ptr[i];
				rowstart = MGS[Cnum].row_ptr[i];
				for (j = 0; j < rownum; ++j) {
					btemp[0] = CBG[Cnum][MGS[Cnum].col_index[rowstart + j] * 3 + 1];
					btemp[1] = CBG[Cnum][MGS[Cnum].col_index[rowstart + j] * 3 + 2];
					btemp[2] = CBG[Cnum][MGS[Cnum].col_index[rowstart + j] * 3 + 3];
					Bvector[ID][i * 3 + 1] += MGS[Cnum].value[rowstart * 6 + j * 6] * btemp[0] + MGS[Cnum].value[rowstart * 6 + j * 6 + 1] * btemp[1] + MGS[Cnum].value[rowstart * 6 + j * 6 + 2] * btemp[2];
					Bvector[ID][i * 3 + 2] += MGS[Cnum].value[rowstart * 6 + j * 6 + 1] * btemp[0] + MGS[Cnum].value[rowstart * 6 + j * 6 + 3] * btemp[1] + MGS[Cnum].value[rowstart * 6 + j * 6 + 4] * btemp[2];
					Bvector[ID][i * 3 + 3] += MGS[Cnum].value[rowstart * 6 + j * 6 + 2] * btemp[0] + MGS[Cnum].value[rowstart * 6 + j * 6 + 4] * btemp[1] + MGS[Cnum].value[rowstart * 6 + j * 6 + 5] * btemp[2];
				}
			}
		}
		

	}


	return nullptr;
}

void* ConvolBvector(Wvector * CBT, Wvector * CBG, BoundaryValue * bd, long StepNum, long MaxN)//串行
{
	long Cnum;
	for (Cnum = 1; Cnum <= StepNum; ++Cnum)
	{
		if (Cnum <= MaxN)
		{
			// [G]
			CBG[Cnum] = 0.0;
			if ((StepNum - Cnum + 1) > 0 && (StepNum - Cnum + 1) < StepNum)
				CBG[Cnum] += bd[StepNum - Cnum + 1].lt;
			if ((StepNum - Cnum) > 0 && (StepNum - Cnum) < StepNum)
				CBG[Cnum].plus(bd[StepNum - Cnum].lt, 2.0);
			if ((StepNum - Cnum - 1) > 0 && (StepNum - Cnum - 1) < StepNum)
				CBG[Cnum] += bd[StepNum - Cnum - 1].lt;
			//T
			CBT[Cnum] = 0.0;
			if ((StepNum - Cnum + 1) > 0 && (StepNum - Cnum + 1) < StepNum)
				CBT[Cnum] += bd[StepNum - Cnum + 1].lu;
			if ((StepNum -Cnum) > 0 && (StepNum - Cnum) < StepNum)
				CBT[Cnum].plus(bd[StepNum - Cnum].lu, 2.0);
			if ((StepNum - Cnum - 1) > 0 && (StepNum - Cnum - 1) < StepNum)
				CBT[Cnum] += bd[StepNum - Cnum - 1].lu;

		}
		else {
			break;
		}
	}

	return nullptr;
}

void* PthConvolBvector(void* arg)//并行
{


	long Cnum;
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = (*PthTemp).targs;
	long TMaxN = (*PthTemp).MaxN;
	long i, j, k, l, pos;
	int Flag;
	long size = (3 * NodeNum) / thread_num;
	long rest = (3 * NodeNum) % thread_num;
	long count = ID * size;
	long end;

	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (Cnum = 1; Cnum <= StepNum; ++Cnum)
	{
		if (Cnum <= TMaxN)
		{

			for (i = count; i < end; ++i)
			{
				CBG[Cnum][i + 1] = 0.0;
			}

			// [G]
			if ((StepNum - Cnum + 1) > 0 && (StepNum - Cnum + 1) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBG[Cnum][i + 1] += bd[StepNum - Cnum + 1].lt[i + 1];
				}
			}
			if ((StepNum - Cnum) > 0 && (StepNum - Cnum) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBG[Cnum][i + 1] += bd[StepNum - Cnum].lt[i + 1] * 2.0;
				}
			}
			if ((StepNum - Cnum - 1) > 0 && (StepNum - Cnum - 1) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBG[Cnum][i + 1] += bd[StepNum - Cnum - 1].lt[i + 1];
				}
			}

			for (i = count; i < end; ++i)
			{
				CBT[Cnum][i + 1] = 0.0;
			}
			//T

			if ((StepNum - Cnum + 1) > 0 && (StepNum - Cnum + 1) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBT[Cnum][i + 1] += bd[StepNum - Cnum + 1].lu[i + 1];
				}
			}
			if ((StepNum - Cnum) > 0 && (StepNum - Cnum) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBT[Cnum][i + 1] += bd[StepNum - Cnum].lu[i + 1] * 2.0;
				}
			}

			if ((StepNum - Cnum - 1) > 0 && (StepNum - Cnum - 1) < StepNum) {
				for (i = count; i < end; ++i)
				{
					CBT[Cnum][i + 1] += bd[StepNum - Cnum - 1].lu[i + 1];
				}
			}

		}
		else {
			break;
		}
	}


	return nullptr;
}

void* PthCollecBA(void* arg)//并行
{
	long Cnum;
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = (*PthTemp).targs;
	long TMaxN = (*PthTemp).MaxN;
	long i, j, k, l, pos;
	int Flag;
	long size = (3 * NodeNum) / thread_num;
	long rest = (3 * NodeNum) % thread_num;
	long count = ID * size;
	long end;
	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)
	{
		for (j = 0; j < thread_num; ++j) {

			Bsum[i + 1] += Bvector[j][i + 1];
			Bvector[j][i + 1] = 0;
	}
	
	}

	return nullptr;
}

void* PthCollecBS(void* arg)
{
	long Cnum;
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = (*PthTemp).targs;
	long TMaxN = (*PthTemp).MaxN;
	long i, j, k, l, pos;
	int Flag;
	long size = (3 * NodeNum) / thread_num;
	long rest = (3 * NodeNum) % thread_num;
	long count = ID * size;
	long end;
	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)
	{
		for (j = 0; j < thread_num; ++j) {

			Bsum[i + 1] -= Bvector[j][i + 1];
			//Bvector[j][i + 1] = 0;
		}

	}

	return nullptr;
}

void* PthBvectorInitial(void* arg)
{
	long Cnum;
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = (*PthTemp).targs;
	long TMaxN = (*PthTemp).MaxN;
	long i, j, k, l, pos;
	int Flag;
	long size = (3 * NodeNum) / thread_num;
	long rest = (3 * NodeNum) % thread_num;
	long count = ID * size;
	long end;
	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)
	{
		for (j = 0; j < thread_num; ++j) {

			Bvector[j][i + 1] = 0;
		}

	}

	return nullptr;
}

void* PthDSEPhyInit(void* arg)
{
	long Cnum;
	PthConvol* PthTemp = (PthConvol*)arg;
	long ID = (*PthTemp).targs;
	long TMaxN = (*PthTemp).MaxN;
	long i, j, k, l, pos;
	int Flag;
	long size = (3 * NodeNum) / thread_num;
	long rest = (3 * NodeNum) % thread_num;
	long count = ID * size;
	long end;
	if (ID < rest) {
		count += ID;
		end = count + size + 1;
	}
	else {
		count += rest;
		end = count + size;
	}

	for (i = count; i < end; ++i)
	{
		for (j = 0; j < thread_num; ++j) {

			Bvector[j][i + 1] = 0;
		}

	}

	return nullptr;
}

