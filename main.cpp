#include "windows.h" // ����������� ��������� ��� ������ � Windows API
#include <vector> // ����������� ��������� ��� ������������� ���������� vector
#include <cmath> // ����������� ��������� ��� �������������� �������
#include <iostream> // ����������� ��������� ��� �����-������
#include <iomanip> // ����������� ��������� ��� ���������� �������� �����-������
#include <locale> // ����������� ��������� ��� ������ � ��������
#include <codecvt> // ����������� ��������� ��� ����������� ���������
#include <sstream> // ����������� ��������� ��� ������ � �������� �����
#include <fstream>

using namespace std; // ������������� ������������ ������������ ����

// ����� ��� ������ � ������������� ��������� � ���������� �� ����������� ��������
class SymQR {
private:
    int dimension; // ����������� �������
    vector<vector<double>> matrix; // ��������� ������ ��� �������� �������
    vector<double> eigenvalues; // ������ ��� �������� ����������� ��������

public:
    // ����������� ������, ���������������� ������� �������� �����������
    SymQR(int dim) : dimension(dim) {
        matrix.resize(dimension, vector<double>(dimension, 0)); // ������������� ������� ������
        eigenvalues.resize(dimension, 0); // ������������� ������� ����������� �������� ������
    }

    void initializeMatrix(const vector<double>& temp) {
        int size = static_cast<int>(sqrt(temp.size()));
        resize(size); // �������� ����������� ����� ����������� �������
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = temp[i * dimension + j];
            }
        }

        if (!isValidMatrix()) {
            throw runtime_error("Matrix is not symmetric.");
        }
    }

    // ����� ��� ������ ������� �� �����
    void readMatrixFromFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Unable to open file.");
        }

        vector<double> temp;
        double value;
        while (file >> value) {
            temp.push_back(value);
        }
        file.close();

        if (temp.empty()) {
            throw runtime_error("The file is empty.");
        }

        int size = static_cast<int>(sqrt(temp.size()));
        if (size * size != temp.size()) {
            throw runtime_error("Invalid matrix format in file. The number of elements is not a perfect square.");
        }

        initializeMatrix(temp); // ���������� ����� ����� ��� ������������� �������
    }

    // ����� ��� ��������� ����������� �������
    void resize(int newDim) {
        if (newDim > 100) {
            throw runtime_error("Dimension exceeds the maximum limit of 100.");
        }
        dimension = newDim; // ��������� �����������
        matrix.resize(dimension, vector<double>(dimension, 0)); // �������� ������ �������
        eigenvalues.resize(dimension, 0); // �������� ������ ������� ����������� ��������
    }

    // ������� ��� ��������, �������� �� ������ ������
    bool isNumber(const wchar_t* str) {
        wchar_t* end;
        wcstod(str, &end); // ����������� ������ � �����
        return *end == L'\0'; // ���� end ��������� �� ����� ������, ������ ��� �����
    }

    // ����� ��� ������ ������� �� ����� �����
    void readMatrix(HWND hWnd) {
        for (int i = 0; i < dimension; i++) { // �������� �� �������
            for (int j = 0; j < dimension; j++) { // �������� �� ��������
                wchar_t buffer[10]; // ����� ��� �������� ������
                GetDlgItemText(hWnd, 100 + i * dimension + j, buffer, 10); // �������� ����� �� ���� �����

                // ���������, �������� �� ��������� �������� ������
                if (!isNumber(buffer)) {
                    throw runtime_error("Invalid input: Please enter numeric values."); // ��������� ����������, ���� ���� ������������
                }

                matrix[i][j] = _wtof(buffer); // ����������� ������ � ����� � ��������� � �������
            }
        }
        // �������� �� �������������� �������
        if (!isValidMatrix()) {
            throw runtime_error("Matrix is not symmetric."); // ��������� ����������, ���� ������� �� �����������
        }
    }

    // ����� ��� ��������, �������� �� ������� ������������
    bool isValidMatrix() {
        for (int i = 0; i < dimension; i++) { // �������� �� �������
            for (int j = 0; j < dimension; j++) { // �������� �� ��������
                if (fabs(matrix[i][j] - matrix[j][i]) > 1e-10) { // ��������� ��������� ���������
                    return false; // ���� �������� �� �����, ���������� false
                }
            }
        }
        return true; // ���� ��� �������� �����, ���������� true
    }

    // ����� ��� QR-���������� �������
    void qrDecomposition(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
        int n = A.size(); // �������� ����������� ������� A
        Q = vector<vector<double>>(n, vector<double>(n, 0)); // ������������� ������� Q
        R = vector<vector<double>>(n, vector<double>(n, 0)); // ������������� ������� R

        for (int j = 0; j < n; j++) { // �������� �� ��������
            vector<double> v(n); // ������ ��� �������� �������� �������
            for (int i = 0; i < n; i++) {
                v[i] = A[i][j]; // �������� ������� �� A � v
            }

            // �����-������� ������� ��� ���������������
            for (int i = 0; i < j; i++) { // �������� �� ���������� ��������
                R[i][j] = 0.0; // ������������� �������� R
                for (int k = 0; k < n; k++) {
                    R[i][j] += Q[k][i] * A[k][j]; // ��������� R
                }
                for (int k = 0; k < n; k++) {
                    v[k] -= R[i][j] * Q[k][i]; // ��������� v
                }
            }

            // ������������ ������� v
            R[j][j] = 0.0; // ������������� ������������� �������� R
            for (int i = 0; i < n; i++) {
                R[j][j] += v[i] * v[i]; // ��������� ������� �����
            }
            R[j][j] = sqrt(R[j][j]); // ����� ������ ��� ��������� �����

            // ��������� ������� Q
            if (R[j][j] != 0) { // ���������, �� ����� �� ������� ����
                for (int i = 0; i < n; i++) {
                    Q[i][j] = v[i] / R[j][j]; // ��������� ������ � ��������� Q
                }
            }
        }
    }

    // ����� ��� ��������� ���� ������
    vector<vector<double>> multiply(const vector<vector<double>>& A, const vector<vector<double>>& B) {
        int n = A.size(); // �������� ����������� ������� A
        vector<vector<double>> C(n, vector<double>(n, 0)); // ������������� �������������� ������� C
        for (int i = 0; i < n; i++) { // �������� �� ������� A
            for (int j = 0; j < n; j++) { // �������� �� �������� B
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j]; // ��������� ������������
                }
            }
        }
        return C; // ���������� ���������
    }

    // ����� ��� ���������� ����������� ��������
    void computeEigenvalues(int maxIterations = 1000, double tolerance = 1e-10) {
        int n = matrix.size(); // �������� ����������� �������
        vector<vector<double>> A = matrix; // �������� �������
        vector<vector<double>> Q(n, vector<double>(n, 0)); // ������������� ������� Q
        vector<vector<double>> R(n, vector<double>(n, 0)); // ������������� ������� R

        for (int k = 0; k < maxIterations; k++) { // ���� ��� ������������� ���������� ��������
            qrDecomposition(A, Q, R); // ��������� QR-����������
            A = multiply(R, Q); // ��������� A

            // �������� �� ����������
            double offDiagonalSum = 0; // ���������� ��� ����� ��������� ��� ������� ���������
            for (int i = 0; i < n; i++) { // �������� �� �������
                for (int j = 0; j < n; j++) { // �������� �� ��������
                    if (i != j) {
                        offDiagonalSum += fabs(A[i][j]); // ��������� �������� ��� ������� ���������
                    }
                }
            }
            if (offDiagonalSum < tolerance) { // ���������, ������ �� ����� ������
                break; // ���� ��, ������� �� �����
            }
        }

        for (int i = 0; i < n; i++) {
            eigenvalues[i] = A[i][i]; // ���������� ����������� �������� �� ������� ���������
        }
    }

    // ����� ��� ��������� ����������� ��������
    vector<double> getEigenvalues() {
        return eigenvalues; // ���������� ������ ����������� ��������
    }

    // ����� ��� �������������� � ������ ����������� ��������
    string displayEigenvalues() {
        ostringstream oss; // ������� ����� ��� ������������ ������
        oss << "Eigenvalues:  \n"; // ��������� ���������
        for (int i = 0; i < dimension; i++) {
            oss << fixed << setprecision(10) << eigenvalues[i] << "; "; // ����������� �����
        }
        return oss.str(); // ���������� ������ � ������������ ����������
    }
};


// ���������� ��������� �������� ����
LRESULT CALLBACK SoftwareMainProcedure(HWND hWnd, UINT msg, WPARAM wp, LPARAM lp);

// ������� ��� ����������� ������ ����
WNDCLASS SymQRClass(HBRUSH BGColor, HCURSOR Cursor, HINSTANCE hInst, HICON Icon, LPCWSTR Name, WNDPROC Procedure) {
    WNDCLASS NWC = { 0 }; // ������������� ��������� WNDCLASS

    NWC.hCursor = Cursor; // ��������� �������
    NWC.hIcon = Icon; // ��������� ������
    NWC.hInstance = hInst; // ��������� ����������
    NWC.lpszClassName = Name; // ��������� ����� ������
    NWC.hbrBackground = BGColor; // ��������� ����
    NWC.lpfnWndProc = Procedure; // ��������� ������� ��������� ���������

    return NWC; // ���������� ������������������ �����
}

// ������� ��� �������� ����� ����� ��� �������
void CreateMatrixInputFields(HWND hWnd, SymQR* qr) {
    int dimension = qr->getEigenvalues().size(); // �������� ������� �����������
    int displayDimension = min(dimension, 10); // ������������ ������������ ����������� �� 10
    int xOffset = 20, yOffset = 40; // ������ �������� ��� ���������������� �����

    // ������� ���� ����� ��� �������
    for (int i = 0; i < displayDimension; i++) { // �������� �� �������
        for (int j = 0; j < displayDimension; j++) { // �������� �� ��������
            CreateWindowW(L"EDIT", L"", // ������� ���� �����
                WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL, // ������� ES_NUMBER, ����� ��������� ���� ����� �����
                xOffset + j * 35, yOffset + i * 35, // �������
                30, 30, // �������
                hWnd, (HMENU)(100 + i * displayDimension + j), NULL, NULL); // ���������� �������������
        }
    }

    // ���������, ���������� �� ������ "���������"
    if (!GetDlgItem(hWnd, 1)) { // ���� ������ �� �������
        // ������� ������ "���������", ���� ��� �� ����������
        CreateWindowW(L"BUTTON", L"���������", // ������� ������
            WS_VISIBLE | WS_CHILD, // ����� ����
            xOffset, yOffset + displayDimension * 30 + 60, // �������
            200, 40, // �������
            hWnd, (HMENU)(1), NULL, NULL); // ���������� �������������
        // �������� ������ "��������� �� �����"
        CreateWindowW(L"BUTTON", L"��������� �� �����", // ������� ������
            WS_VISIBLE | WS_CHILD, // ����� ����
            xOffset + 200, yOffset + displayDimension * 30 + 60, // �������
            200, 40, // �������
            hWnd, (HMENU)(2), NULL, NULL); // ���������� ������������� (2)
    }
}

// ������� ��� ���������� ����� ����� �������
void UpdateMatrixInputFields(HWND hWnd, SymQR* qr) {
    // ������� ������ ���� �����
    for (int i = 0; i < 10; i++) { // ������������ ����������� 10
        for (int j = 0; j < 10; j++) {
            HWND hEdit = GetDlgItem(hWnd, 100 + i * 10 + j); // �������� ��������� �� ���� �����
            if (hEdit) {
                DestroyWindow(hEdit); // ���������� ������ ���� �����
            }
        }
    }

    // ������� ����� ���� ����� � ����� ������������
    CreateMatrixInputFields(hWnd, qr); // ������� ����� ���� �����
}

// ������� ��� ����������� ����������� ��������
void ShowEigenvalues(HWND hWnd, SymQR* qr, HWND hResultTextBox) {
    string result = qr->displayEigenvalues(); // �������� ������ � ������������ ����������
    // ����������� ����������� �������� � ��������� ����
    SetWindowTextA(hResultTextBox, result.c_str()); // ������������� ����� � ��������� ����
}

// �������� ������� ����������
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR args, int ncmdshow) {
    // ����������� ������ ����
    WNDCLASS SoftwareMainClass = SymQRClass((HBRUSH)COLOR_WINDOW, LoadCursor(NULL, IDC_ARROW), hInst, LoadIcon(NULL, IDI_QUESTION),
        L"SymQR", SoftwareMainProcedure); // �������� ������ ����

    if (!RegisterClassW(&SoftwareMainClass)) { return -1; } // ���� ����������� �� �������, �������
    MSG SoftwareMainMessage = { 0 }; // ������������� ��������� MSG

    // �������� �������� ����
    CreateWindow(L"SymQR", L"QR Matrix", WS_OVERLAPPEDWINDOW | WS_VISIBLE, 10, 10, 800, 800, NULL, NULL, hInst, NULL);

    // ������� ���� ��������� ���������
    while (GetMessage(&SoftwareMainMessage, NULL, NULL, NULL)) { // �������� ���������
        TranslateMessage(&SoftwareMainMessage); // ��������� ���������
        DispatchMessageW(&SoftwareMainMessage); // ������������ ���������
    }
    return 0; // ���������� 0 ��� ����������
}

// ������� ��������� ��������� ��� �������� ����
LRESULT CALLBACK SoftwareMainProcedure(HWND hWnd, UINT msg, WPARAM wp, LPARAM lp) {
    static SymQR* qr = nullptr; // ��������� �� ������ SymQR
    static HWND hResultTextBox = NULL; // ��������� �� ��������� ���� ��� ������ ����������� ��������
    static HWND hDimensionComboBox = NULL; // ��������� �� ���������� ������ ��� ������ �����������
    static HWND hDimensionLabel = NULL; // ��������� �� ����������� ����� ��� �����������
    static HWND hDimensionLabel1 = NULL; // ��������� �� ����������� ����� ��� �������
    int N = 100; // ��������� ����������� �������

    switch (msg) { // ��������� ���������
    case WM_CREATE: // ��������� � �������� ����
        qr = new SymQR(N); // ������� ������ SymQR ��� 10x10 �������
        CreateMatrixInputFields(hWnd, qr); // ������� ���� ����� ��� �������

        // �������� ���������� ���� ��� ����������� ����������� ��������
        hResultTextBox = CreateWindowW(L"EDIT", L"",
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_MULTILINE | ES_READONLY | WS_VSCROLL,
            20, 600, 440, 150, // ������ ���������� � �������
            hWnd, NULL, NULL, NULL);

        // �������� ���������� ���� ��� ����������� ����� �������
        hDimensionLabel1 = CreateWindowW(L"STATIC", L"�������:",
            WS_CHILD | WS_VISIBLE,
            20, 10, 300, 20, // ������ ���������� � �������
            hWnd, NULL, NULL, NULL);

        // �������� ������������ ������ "�������� ����������� �������"
        hDimensionLabel = CreateWindowW(L"STATIC", L"�������� ����������� �������:",
            WS_CHILD | WS_VISIBLE,
            20, 470, 300, 20, // ������ ���������� � �������
            hWnd, NULL, NULL, NULL);

        // �������� ����������� ������ ��� ������ ����������� �������
        hDimensionComboBox = CreateWindowW(L"COMBOBOX", NULL,
            WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST,
            20, 500, 100, 100,
            hWnd, NULL, NULL, NULL);

        // ���������� ��������� � ���������� ������
        for (int i = 2; i <= 10; i++) { // ��������� ����������� �� 2 �� 10
            wchar_t buffer[10]; // ����� ��� �������� ������
            swprintf(buffer, 10, L"%d", i); // ����������� ����� � ������
            SendMessage(hDimensionComboBox, CB_ADDSTRING, 0, (LPARAM)buffer); // ��������� �����������
        }
        SendMessage(hDimensionComboBox, CB_SETCURSEL, 0, 0); // ��������� ���������� ��������
        break;

    case WM_COMMAND: // ��������� ������ (��������, ������� ������)
        if (LOWORD(wp) == 1) { // ID ������ "���������"
            try {
                qr->readMatrix(hWnd); // ������ ������� �� ����� �����
                qr->computeEigenvalues(); // ��������� ����������� ��������
                ShowEigenvalues(hWnd, qr, hResultTextBox); // ���������� ����������� ��������
            }
            catch (const runtime_error& e) { // ��������� ����������
                // �������������� char* � wchar_t*
                wstring_convert<codecvt_utf8<wchar_t>, wchar_t> converter; // ��������� ��� ��������������
                wstring wideError = converter.from_bytes(e.what()); // ����������� ��������� �� ������

                MessageBox(hWnd, wideError.c_str(), L"Error", MB_OK | MB_ICONERROR); // ���������� ��������� �� ������
            }
        }
        else if (LOWORD(wp) == 2) { // ID ������ "��������� �� �����"
            try {
                qr->readMatrixFromFile("file.txt"); // ������ ������� �� �����
                qr->computeEigenvalues(); // ��������� ����������� ��������
                ShowEigenvalues(hWnd, qr, hResultTextBox); // ���������� ����������� ��������

                // ��������� ���� ����� � ����� ������������
                UpdateMatrixInputFields(hWnd, qr); // ��������� ���� �����
            }
            catch (const runtime_error& e) {
                wstring_convert<codecvt_utf8<wchar_t>, wchar_t> converter;
                wstring wideError = converter.from_bytes(e.what());
                MessageBox(hWnd, wideError.c_str(), L"Error", MB_OK | MB_ICONERROR);
            }
        }
        else if (HIWORD(wp) == CBN_SELCHANGE && (HWND)lp == hDimensionComboBox) { // ���� ��������� ����� � ����������
            int selectedDimension = SendMessage(hDimensionComboBox, CB_GETCURSEL, 0, 0) + 2; // �������� ��������� �����������
            qr->resize(selectedDimension); // �������� ����������� �������
            UpdateMatrixInputFields(hWnd, qr); // ��������� ���� �����
        }
        break;

    case WM_DESTROY: // ��������� � �������� ����
        delete qr; // ����������� ������
        PostQuitMessage(0); // ��������� ����������
        break;
    default:
        return DefWindowProc(hWnd, msg, wp, lp); // ��������� ����������� ���������
    }
    return 0; // ���������� 0
}