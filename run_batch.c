/* Windows C batch runner for DBEM1. Build: cl /nologo /W4 /O2 run_batch.c */
#define _CRT_SECURE_NO_WARNINGS
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PATH_SIZE 2048
#define NAME_SIZE 512

typedef struct {
    char requested[NAME_SIZE];
    char relative[NAME_SIZE];
    char card[PATH_SIZE];
    char error[NAME_SIZE];
} Case;

static volatile LONG stopped = 0;

static BOOL WINAPI on_ctrl(DWORD event)
{
    if (event == CTRL_C_EVENT || event == CTRL_BREAK_EVENT || event == CTRL_CLOSE_EVENT) {
        InterlockedExchange(&stopped, 1);
        return TRUE;
    }
    return FALSE;
}

static int exists(const char *path)
{
    DWORD attr = GetFileAttributesA(path);
    return attr != INVALID_FILE_ATTRIBUTES && !(attr & FILE_ATTRIBUTE_DIRECTORY);
}

static int is_directory(const char *path)
{
    DWORD attr = GetFileAttributesA(path);
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY);
}

static int find_project_root(char *root, size_t capacity)
{
    char *slash;
    DWORD length = GetModuleFileNameA(NULL, root, (DWORD)capacity);
    if (!length || length >= capacity) return 0;
    slash = strrchr(root, '\\');
    if (!slash) return 0;
    *slash = 0;
    for (;;) {
        if (!SetCurrentDirectoryA(root)) return 0;
        if (is_directory("DBEM1\\input") && exists("DBEM1\\BEM_DATACARD.DAT")) return 1;
        slash = strrchr(root, '\\');
        if (!slash || slash == root + 2) break;
        *slash = 0;
    }
    return 0;
}

static int finish(int status, int own_console)
{
    if (own_console) {
        printf("\nPress Enter to close this window...");
        fflush(stdout);
        getchar();
    }
    return status;
}

static int missing_dependency(const char *exe, char *missing, size_t capacity)
{
    static const char *dlls[] = {
        "pthreadVC2.dll", "MSVCP140.dll", "VCOMP140.DLL",
        "VCRUNTIME140_1.dll", "VCRUNTIME140.dll"
    };
    char full[PATH_SIZE], local[PATH_SIZE], *slash;
    size_t i;
    DWORD length = GetFullPathNameA(exe, sizeof(full), full, NULL);
    if (!length || length >= sizeof(full)) return 0;
    slash = strrchr(full, '\\');
    if (!slash) return 0;
    *slash = 0;
    for (i = 0; i < sizeof(dlls) / sizeof(dlls[0]); ++i) {
        if (snprintf(local, sizeof(local), "%s\\%s", full, dlls[i]) >= (int)sizeof(local)) continue;
        if (!exists(local) && !SearchPathA(NULL, dlls[i], NULL, 0, NULL, NULL)) {
            strcpy_s(missing, capacity, dlls[i]);
            return 1;
        }
    }
    return 0;
}

static void trim(char *text)
{
    size_t n;
    char *start = text;
    while (*start && isspace((unsigned char)*start)) start++;
    if (start != text) memmove(text, start, strlen(start) + 1);
    n = strlen(text);
    while (n && isspace((unsigned char)text[n - 1])) text[--n] = 0;
}

static int valid_name(char *name)
{
    char *part = name;
    size_t i, n = strlen(name);
    if (!n || n > 48) return 0; /* fscanf_s(..., 49) needs one byte for NUL. */
    for (i = 0; i <= n; ++i) {
        unsigned char ch = (unsigned char)name[i];
        if (ch == '/' || ch == '\\' || ch == 0) {
            if (name + i == part || (name + i - part == 1 && part[0] == '.') ||
                (name + i - part == 2 && part[0] == '.' && part[1] == '.')) return 0;
            part = name + i + 1;
            if (ch == '/') name[i] = '\\';
        } else if (ch < 33 || ch >= 127 || strchr("<>:\"|?*", ch)) {
            return 0;
        }
    }
    return 1;
}

static void find_models(const char *dir, const char *filename, char *first, int *count)
{
    WIN32_FIND_DATAA entry;
    HANDLE search;
    char pattern[PATH_SIZE], path[PATH_SIZE];
    if (snprintf(pattern, sizeof(pattern), "%s\\*", dir) >= (int)sizeof(pattern)) return;
    search = FindFirstFileA(pattern, &entry);
    if (search == INVALID_HANDLE_VALUE) return;
    do {
        if (!strcmp(entry.cFileName, ".") || !strcmp(entry.cFileName, "..")) continue;
        if (snprintf(path, sizeof(path), "%s\\%s", dir, entry.cFileName) >= (int)sizeof(path)) continue;
        if (entry.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            if (!(entry.dwFileAttributes & FILE_ATTRIBUTE_REPARSE_POINT))
                find_models(path, filename, first, count);
        } else if (!_stricmp(entry.cFileName, filename)) {
            if (++*count == 1) strcpy_s(first, PATH_SIZE, path);
        }
    } while (FindNextFileA(search, &entry));
    FindClose(search);
}

static void resolve_case(Case *item)
{
    char name[NAME_SIZE], model[PATH_SIZE], path[PATH_SIZE];
    int count = 0;
    strcpy_s(name, sizeof(name), item->requested);
    if (!valid_name(name)) {
        strcpy_s(item->error, sizeof(item->error), "invalid name (relative ASCII path, max 48 bytes, no spaces)");
        return;
    }
    if (strrchr(name, '\\')) {
        if (snprintf(model, sizeof(model), "DBEM1\\input\\%s.model", name) >= (int)sizeof(model)) {
            strcpy_s(item->error, sizeof(item->error), "model path too long");
            return;
        }
        count = exists(model) ? 1 : 0;
    } else {
        if (snprintf(path, sizeof(path), "%s.model", name) >= (int)sizeof(path)) return;
        find_models("DBEM1\\input", path, model, &count);
    }
    if (count > 1) {
        strcpy_s(item->error, sizeof(item->error), "duplicate model; specify its path below DBEM1/input");
        return;
    }
    if (!count) {
        strcpy_s(item->error, sizeof(item->error), "missing .model");
        return;
    }
    model[strlen(model) - 6] = 0; /* Remove .model. */
    if (snprintf(path, sizeof(path), "%s.bd", model) >= (int)sizeof(path) || !exists(path)) {
        strcpy_s(item->error, sizeof(item->error), "missing .bd");
        return;
    }
    if (snprintf(item->card, sizeof(item->card), "%s.DATACARD", model) >= (int)sizeof(item->card)
        || !exists(item->card)) {
        strcpy_s(item->error, sizeof(item->error), "missing .DATACARD");
        return;
    }
    /* model starts with DBEM1\input\; DBEM1 then prepends .\input\ itself. */
    strcpy_s(item->relative, sizeof(item->relative), model + strlen("DBEM1\\input\\"));
}

static int load_cases(const char *list, Case **items, size_t *count)
{
    FILE *in = fopen(list, "rb");
    Case *cases = NULL;
    size_t used = 0, capacity = 0;
    char line[NAME_SIZE];
    if (!in) return 0;
    while (fgets(line, sizeof(line), in)) {
        Case *grown;
        if (!strchr(line, '\n') && !feof(in)) {
            fprintf(stderr, "case list line too long\n");
            free(cases); fclose(in); return 0;
        }
        if (!used && (unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB
            && (unsigned char)line[2] == 0xBF) memmove(line, line + 3, strlen(line + 3) + 1);
        trim(line);
        if (!line[0] || line[0] == '#') continue;
        if (used == capacity) {
            capacity = capacity ? capacity * 2 : 16;
            grown = (Case *)realloc(cases, capacity * sizeof(Case));
            if (!grown) { free(cases); fclose(in); return 0; }
            cases = grown;
        }
        memset(&cases[used], 0, sizeof(Case));
        strcpy_s(cases[used].requested, sizeof(cases[used].requested), line);
        resolve_case(&cases[used]);
        used++;
    }
    fclose(in);
    *items = cases;
    *count = used;
    return used != 0;
}

static int write_card(const char *source, const char *name)
{
    FILE *in = NULL, *out = NULL;
    char *data = NULL;
    long size;
    size_t begin = 0, end;
    int ok = 0;
    in = fopen(source, "rb");
    if (!in || fseek(in, 0, SEEK_END) || (size = ftell(in)) < 0 || fseek(in, 0, SEEK_SET)) goto done;
    data = (char *)malloc((size_t)size + 1);
    if (!data || fread(data, 1, (size_t)size, in) != (size_t)size) goto done;
    while (begin < (size_t)size && isspace((unsigned char)data[begin])) begin++;
    if (begin == (size_t)size || (size >= 3 && !memcmp(data, "\xef\xbb\xbf", 3))) goto done;
    end = begin;
    while (end < (size_t)size && !isspace((unsigned char)data[end])) end++;
    out = fopen("DBEM1\\BEM_DATACARD.DAT.batch.tmp", "wb");
    if (!out || fwrite(data, 1, begin, out) != begin ||
        fwrite(name, 1, strlen(name), out) != strlen(name) ||
        fwrite(data + end, 1, (size_t)size - end, out) != (size_t)size - end) goto done;
    if (fclose(out)) { out = NULL; goto done; }
    out = NULL;
    ok = MoveFileExA("DBEM1\\BEM_DATACARD.DAT.batch.tmp", "DBEM1\\BEM_DATACARD.DAT",
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) != 0;
done:
    if (out) fclose(out);
    if (in) fclose(in);
    free(data);
    DeleteFileA("DBEM1\\BEM_DATACARD.DAT.batch.tmp");
    return ok;
}

static void csv_field(FILE *out, const char *value)
{
    const char *p;
    fputc('"', out);
    for (p = value; *p; ++p) { if (*p == '"') fputc('"', out); fputc(*p, out); }
    fputc('"', out);
}

static int changed_log(const FILETIME *before)
{
    WIN32_FILE_ATTRIBUTE_DATA data;
    return GetFileAttributesExA("output\\log.txt", GetFileExInfoStandard, &data)
        && (!before || CompareFileTime(before, &data.ftLastWriteTime) != 0);
}

static int launch_case(const char *exe, DWORD *code, DWORD *error)
{
    char full_exe[PATH_SIZE], full_dir[PATH_SIZE], command[PATH_SIZE + 4];
    STARTUPINFOA startup;
    PROCESS_INFORMATION process;
    DWORD wait_result;
    *error = 0;
    if (!GetFullPathNameA(exe, sizeof(full_exe), full_exe, NULL) ||
        !GetFullPathNameA("DBEM1", sizeof(full_dir), full_dir, NULL)) { *error = GetLastError(); return 0; }
    if (snprintf(command, sizeof(command), "\"%s\"", full_exe) >= (int)sizeof(command)) return 0;
    ZeroMemory(&startup, sizeof(startup));
    ZeroMemory(&process, sizeof(process));
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    startup.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    startup.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    fflush(stdout);
    fflush(stderr);
    if (!CreateProcessA(full_exe, command, NULL, NULL, TRUE, 0, NULL, full_dir, &startup, &process)) {
        *error = GetLastError();
        return 0;
    }
    CloseHandle(process.hThread);
    for (;;) {
        wait_result = WaitForSingleObject(process.hProcess, 200);
        if (wait_result == WAIT_OBJECT_0) break;
        if (wait_result != WAIT_TIMEOUT || InterlockedCompareExchange(&stopped, 0, 0)) {
            TerminateProcess(process.hProcess, 130);
            WaitForSingleObject(process.hProcess, INFINITE);
            break;
        }
    }
    if (!GetExitCodeProcess(process.hProcess, code)) *code = (DWORD)-1;
    CloseHandle(process.hProcess);
    return 1;
}

int main(int argc, char **argv)
{
    char root[PATH_SIZE], list[PATH_SIZE] = "cases.txt", exe[PATH_SIZE] = "x64\\Release\\DBEM1.exe";
    char batch[PATH_SIZE], backup[PATH_SIZE], lock_text[PATH_SIZE + 80];
    Case *cases = NULL;
    size_t count = 0, i;
    int dry = 0, failed = 0, exit_status = 0, backup_ready = 0, seen_list = 0;
    HANDLE lock = INVALID_HANDLE_VALUE;
    SYSTEMTIME now;
    DWORD written, code, console_ids[2], launch_error;
    int own_console = GetConsoleProcessList(console_ids, 2) == 1;
    char missing_dll[64];
    FILE *summary = NULL;
    if (!find_project_root(root, sizeof(root))) {
        fprintf(stderr, "cannot find project root above run_batch.exe (need DBEM1/input and DBEM1/BEM_DATACARD.DAT)\n");
        return finish(2, own_console);
    }
    for (i = 1; i < (size_t)argc; ++i) {
        if (!strcmp(argv[i], "--dry-run")) dry = 1;
        else if (!strcmp(argv[i], "--exe") && i + 1 < (size_t)argc) strcpy_s(exe, sizeof(exe), argv[++i]);
        else if (argv[i][0] != '-' && !seen_list) {
            strcpy_s(list, sizeof(list), argv[i]); seen_list = 1;
        }
        else { fprintf(stderr, "usage: run_batch.exe [cases.txt] [--exe path] [--dry-run]\n");
               return finish(2, own_console); }
    }
    {
        const char *required[] = { exe, "DBEM1\\gaussposition", "DBEM1\\gaussweight",
                                   "DBEM1\\BEM_DATACARD.DAT" };
        size_t j;
        int missing_file = 0;
        for (j = 0; j < sizeof(required) / sizeof(required[0]); ++j) {
            if (!exists(required[j])) { fprintf(stderr, "missing file: %s\n", required[j]); missing_file = 1; }
        }
        if (missing_file) {
            fprintf(stderr, "project folder: %s\n", root);
            return finish(2, own_console);
        }
    }
    if (!load_cases(list, &cases, &count)) {
        fprintf(stderr, "cannot read cases or list is empty: %s (program folder: %s)\n", list, root);
        return finish(2, own_console);
    }
    for (i = 0; i < count; ++i) {
        printf("[%lu/%lu] %s: %s\n", (unsigned long)(i + 1), (unsigned long)count,
               cases[i].requested, cases[i].error[0] ? cases[i].error : "ready");
        if (cases[i].error[0]) failed++;
    }
    if (!_stricmp(exe, "x64\\Release\\DBEM1.exe") &&
        missing_dependency(exe, missing_dll, sizeof(missing_dll))) {
        fprintf(stderr, "DBEM1.exe needs %s, but Windows cannot find it.\n", missing_dll);
        fprintf(stderr, "For pthreadVC2.dll, place the x64 DLL beside DBEM1.exe; for Visual C++ DLLs, install the x64 runtime.\n");
        free(cases); return finish(2, own_console);
    }
    if (dry) {
        printf("checked: %lu ready, %d failed; no files changed\n", (unsigned long)count - failed, failed);
        free(cases); return finish(failed ? 1 : 0, own_console);
    }
    CreateDirectoryA("output", NULL);
    lock = CreateFileA("output\\.dbem_batch.lock", GENERIC_WRITE, 0, NULL, CREATE_NEW,
                       FILE_ATTRIBUTE_NORMAL, NULL);
    if (lock == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "batch lock exists or cannot be created: output/.dbem_batch.lock\n");
        free(cases); return finish(2, own_console);
    }
    SetConsoleCtrlHandler(on_ctrl, TRUE);
    GetLocalTime(&now);
    snprintf(batch, sizeof(batch), "output\\batch_%04u%02u%02u_%02u%02u%02u_%03u_%lu",
             now.wYear, now.wMonth, now.wDay, now.wHour, now.wMinute, now.wSecond, now.wMilliseconds,
             (unsigned long)GetCurrentProcessId());
    snprintf(backup, sizeof(backup), "%s\\BEM_DATACARD.original.DAT", batch);
    snprintf(lock_text, sizeof(lock_text), "pid=%lu\r\nbackup=%s\r\n",
             (unsigned long)GetCurrentProcessId(), backup);
    WriteFile(lock, lock_text, (DWORD)strlen(lock_text), &written, NULL);
    CloseHandle(lock);
    if (!CreateDirectoryA(batch, NULL) || !CopyFileA("DBEM1\\BEM_DATACARD.DAT", backup, TRUE)) {
        fprintf(stderr, "cannot create batch folder or back up card\n");
        exit_status = 2; goto cleanup;
    }
    backup_ready = 1;
    {
        char summary_path[PATH_SIZE];
        snprintf(summary_path, sizeof(summary_path), "%s\\summary.csv", batch);
        summary = fopen(summary_path, "wb");
    }
    if (!summary) { fprintf(stderr, "cannot create summary.csv\n"); exit_status = 2; goto cleanup; }
    fprintf(summary, "number,requested,matched,status,return_code,message\r\n");
    failed = 0;
    for (i = 0; i < count; ++i) {
        const char *status = "SUCCESS", *message = "";
        char launch_message[NAME_SIZE];
        int started = 0;
        FILETIME before, *before_ptr = NULL;
        WIN32_FILE_ATTRIBUTE_DATA log_info;
        char log_copy[PATH_SIZE];
        code = (DWORD)-1;
        if (InterlockedCompareExchange(&stopped, 0, 0)) break;
        if (cases[i].error[0]) { status = "FAILED"; message = cases[i].error; }
        else {
            if (GetFileAttributesExA("output\\log.txt", GetFileExInfoStandard, &log_info)) {
                before = log_info.ftLastWriteTime; before_ptr = &before;
            }
            if (!write_card(cases[i].card, cases[i].relative)) {
                status = "FAILED"; message = "cannot prepare BEM_DATACARD.DAT";
            } else if (!launch_case(exe, &code, &launch_error)) {
                snprintf(launch_message, sizeof(launch_message),
                         "cannot start DBEM1.exe (Win32 error %lu)", (unsigned long)launch_error);
                status = "FAILED"; message = launch_message;
            } else {
                started = 1;
                if (InterlockedCompareExchange(&stopped, 0, 0)) {
                    status = "INTERRUPTED"; message = "user interrupted";
                } else if (code != 0 && code != 1) {
                    status = "FAILED";
                    if (code == 0xC0000135) message = "DBEM1 missing a DLL (0xC0000135)";
                    else if (code == 0xC000007B) message = "DBEM1 DLL architecture mismatch (0xC000007B)";
                    else {
                        snprintf(launch_message, sizeof(launch_message),
                                 "DBEM1 returned %lu (0x%08lX); see DBEM log",
                                 (unsigned long)code, (unsigned long)code);
                        message = launch_message;
                    }
                }
            }
            if (started && changed_log(before_ptr)) {
                snprintf(log_copy, sizeof(log_copy), "%s\\%04lu_dbem.log", batch, (unsigned long)(i + 1));
                CopyFileA("output\\log.txt", log_copy, FALSE);
            }
        }
        fprintf(summary, "%lu,", (unsigned long)(i + 1));
        csv_field(summary, cases[i].requested); fputc(',', summary);
        csv_field(summary, cases[i].relative); fputc(',', summary);
        csv_field(summary, status); fputc(',', summary);
        if (started) fprintf(summary, "%lu", (unsigned long)code);
        fputc(',', summary); csv_field(summary, message); fputs("\r\n", summary);
        fflush(summary);
        printf("[%lu/%lu] %s: %s%s%s\n", (unsigned long)(i + 1), (unsigned long)count,
               cases[i].requested, status, *message ? " - " : "", message);
        if (strcmp(status, "SUCCESS")) failed++;
    }
    printf("finished: %lu successful, %d failed/interrupted; logs: %s\n",
           (unsigned long)i - failed, failed, batch);
    exit_status = InterlockedCompareExchange(&stopped, 0, 0) ? 130 : (failed ? 1 : 0);
cleanup:
    if (summary) fclose(summary);
    if (backup_ready && !CopyFileA(backup, "DBEM1\\BEM_DATACARD.DAT", FALSE)) {
        fprintf(stderr, "CARD RESTORE FAILED. Restore from %s; lock retained.\n", backup);
        exit_status = 2;
    } else {
        DeleteFileA("output\\.dbem_batch.lock");
    }
    SetConsoleCtrlHandler(on_ctrl, FALSE);
    free(cases);
    return finish(exit_status, own_console);
}
