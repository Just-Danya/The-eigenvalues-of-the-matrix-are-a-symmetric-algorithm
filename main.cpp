#include "windows.h" // Подключение заголовка для работы с Windows API
#include <vector> // Подключение заголовка для использования контейнера vector
#include <cmath> // Подключение заголовка для математических функций
#include <iostream> // Подключение заголовка для ввода-вывода
#include <iomanip> // Подключение заголовка для управления форматом ввода-вывода
#include <locale> // Подключение заголовка для работы с локалями
#include <codecvt> // Подключение заголовка для конвертации кодировок
#include <sstream> // Подключение заголовка для работы с потоками строк
#include <fstream>

using namespace std; // Использование стандартного пространства имен

// Класс для работы с симметричными матрицами и вычисления их собственных значений
class SymQR {
private:
    int dimension; // Размерность матрицы
    vector<vector<double>> matrix; // Двумерный вектор для хранения матрицы
    vector<double> eigenvalues; // Вектор для хранения собственных значений

public:
    // Конструктор класса, инициализирующий матрицу заданной размерности
    SymQR(int dim) : dimension(dim) {
        matrix.resize(dimension, vector<double>(dimension, 0)); // Инициализация матрицы нулями
        eigenvalues.resize(dimension, 0); // Инициализация вектора собственных значений нулями
    }

    void initializeMatrix(const vector<double>& temp) {
        int size = static_cast<int>(sqrt(temp.size()));
        resize(size); // Изменяем размерность перед заполнением матрицы
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = temp[i * dimension + j];
            }
        }

        if (!isValidMatrix()) {
            throw runtime_error("Matrix is not symmetric.");
        }
    }

    // Метод для чтения матрицы из файла
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

        initializeMatrix(temp); // Используем новый метод для инициализации матрицы
    }

    // Метод для изменения размерности матрицы
    void resize(int newDim) {
        if (newDim > 100) {
            throw runtime_error("Dimension exceeds the maximum limit of 100.");
        }
        dimension = newDim; // Обновляем размерность
        matrix.resize(dimension, vector<double>(dimension, 0)); // Изменяем размер матрицы
        eigenvalues.resize(dimension, 0); // Изменяем размер вектора собственных значений
    }

    // Функция для проверки, является ли строка числом
    bool isNumber(const wchar_t* str) {
        wchar_t* end;
        wcstod(str, &end); // Преобразуем строку в число
        return *end == L'\0'; // Если end указывает на конец строки, значит это число
    }

    // Метод для чтения матрицы из полей ввода
    void readMatrix(HWND hWnd) {
        for (int i = 0; i < dimension; i++) { // Проходим по строкам
            for (int j = 0; j < dimension; j++) { // Проходим по столбцам
                wchar_t buffer[10]; // Буфер для хранения строки
                GetDlgItemText(hWnd, 100 + i * dimension + j, buffer, 10); // Получаем текст из поля ввода

                // Проверяем, является ли введенное значение числом
                if (!isNumber(buffer)) {
                    throw runtime_error("Invalid input: Please enter numeric values."); // Генерация исключения, если ввод некорректный
                }

                matrix[i][j] = _wtof(buffer); // Преобразуем строку в число и сохраняем в матрицу
            }
        }
        // Проверка на симметричность матрицы
        if (!isValidMatrix()) {
            throw runtime_error("Matrix is not symmetric."); // Генерация исключения, если матрица не симметрична
        }
    }

    // Метод для проверки, является ли матрица симметричной
    bool isValidMatrix() {
        for (int i = 0; i < dimension; i++) { // Проходим по строкам
            for (int j = 0; j < dimension; j++) { // Проходим по столбцам
                if (fabs(matrix[i][j] - matrix[j][i]) > 1e-10) { // Проверяем равенство элементов
                    return false; // Если элементы не равны, возвращаем false
                }
            }
        }
        return true; // Если все элементы равны, возвращаем true
    }

    // Метод для QR-разложения матрицы
    void qrDecomposition(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
        int n = A.size(); // Получаем размерность матрицы A
        Q = vector<vector<double>>(n, vector<double>(n, 0)); // Инициализация матрицы Q
        R = vector<vector<double>>(n, vector<double>(n, 0)); // Инициализация матрицы R

        for (int j = 0; j < n; j++) { // Проходим по столбцам
            vector<double> v(n); // Вектор для хранения текущего столбца
            for (int i = 0; i < n; i++) {
                v[i] = A[i][j]; // Копируем столбец из A в v
            }

            // Грамм-Шмидтов процесс для ортогонализации
            for (int i = 0; i < j; i++) { // Проходим по предыдущим столбцам
                R[i][j] = 0.0; // Инициализация элемента R
                for (int k = 0; k < n; k++) {
                    R[i][j] += Q[k][i] * A[k][j]; // Вычисляем R
                }
                for (int k = 0; k < n; k++) {
                    v[k] -= R[i][j] * Q[k][i]; // Обновляем v
                }
            }

            // Нормализация вектора v
            R[j][j] = 0.0; // Инициализация диагонального элемента R
            for (int i = 0; i < n; i++) {
                R[j][j] += v[i] * v[i]; // Вычисляем квадрат нормы
            }
            R[j][j] = sqrt(R[j][j]); // Берем корень для получения нормы

            // Заполняем матрицу Q
            if (R[j][j] != 0) { // Проверяем, не равен ли элемент нулю
                for (int i = 0; i < n; i++) {
                    Q[i][j] = v[i] / R[j][j]; // Нормируем вектор и заполняем Q
                }
            }
        }
    }

    // Метод для умножения двух матриц
    vector<vector<double>> multiply(const vector<vector<double>>& A, const vector<vector<double>>& B) {
        int n = A.size(); // Получаем размерность матрицы A
        vector<vector<double>> C(n, vector<double>(n, 0)); // Инициализация результирующей матрицы C
        for (int i = 0; i < n; i++) { // Проходим по строкам A
            for (int j = 0; j < n; j++) { // Проходим по столбцам B
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j]; // Суммируем произведения
                }
            }
        }
        return C; // Возвращаем результат
    }

    // Метод для вычисления собственных значений
    void computeEigenvalues(int maxIterations = 1000, double tolerance = 1e-10) {
        int n = matrix.size(); // Получаем размерность матрицы
        vector<vector<double>> A = matrix; // Копируем матрицу
        vector<vector<double>> Q(n, vector<double>(n, 0)); // Инициализация матрицы Q
        vector<vector<double>> R(n, vector<double>(n, 0)); // Инициализация матрицы R

        for (int k = 0; k < maxIterations; k++) { // Цикл для максимального количества итераций
            qrDecomposition(A, Q, R); // Выполняем QR-разложение
            A = multiply(R, Q); // Обновляем A

            // Проверка на сходимость
            double offDiagonalSum = 0; // Переменная для суммы элементов вне главной диагонали
            for (int i = 0; i < n; i++) { // Проходим по строкам
                for (int j = 0; j < n; j++) { // Проходим по столбцам
                    if (i != j) {
                        offDiagonalSum += fabs(A[i][j]); // Суммируем элементы вне главной диагонали
                    }
                }
            }
            if (offDiagonalSum < tolerance) { // Проверяем, меньше ли сумма порога
                break; // Если да, выходим из цикла
            }
        }

        for (int i = 0; i < n; i++) {
            eigenvalues[i] = A[i][i]; // Записываем собственные значения из главной диагонали
        }
    }

    // Метод для получения собственных значений
    vector<double> getEigenvalues() {
        return eigenvalues; // Возвращаем вектор собственных значений
    }

    // Метод для форматирования и вывода собственных значений
    string displayEigenvalues() {
        ostringstream oss; // Создаем поток для формирования строки
        oss << "Eigenvalues:  \n"; // Добавляем заголовок
        for (int i = 0; i < dimension; i++) {
            oss << fixed << setprecision(10) << eigenvalues[i] << "; "; // Форматируем вывод
        }
        return oss.str(); // Возвращаем строку с собственными значениями
    }
};


// Обработчик сообщений главного окна
LRESULT CALLBACK SoftwareMainProcedure(HWND hWnd, UINT msg, WPARAM wp, LPARAM lp);

// Функция для регистрации класса окна
WNDCLASS SymQRClass(HBRUSH BGColor, HCURSOR Cursor, HINSTANCE hInst, HICON Icon, LPCWSTR Name, WNDPROC Procedure) {
    WNDCLASS NWC = { 0 }; // Инициализация структуры WNDCLASS

    NWC.hCursor = Cursor; // Установка курсора
    NWC.hIcon = Icon; // Установка иконки
    NWC.hInstance = hInst; // Установка экземпляра
    NWC.lpszClassName = Name; // Установка имени класса
    NWC.hbrBackground = BGColor; // Установка фона
    NWC.lpfnWndProc = Procedure; // Установка функции обработки сообщений

    return NWC; // Возвращаем зарегистрированный класс
}

// Функция для создания полей ввода для матрицы
void CreateMatrixInputFields(HWND hWnd, SymQR* qr) {
    int dimension = qr->getEigenvalues().size(); // Получаем текущую размерность
    int displayDimension = min(dimension, 10); // Ограничиваем отображаемую размерность до 10
    int xOffset = 20, yOffset = 40; // Задаем смещения для позиционирования полей

    // Создаем поля ввода для матрицы
    for (int i = 0; i < displayDimension; i++) { // Проходим по строкам
        for (int j = 0; j < displayDimension; j++) { // Проходим по столбцам
            CreateWindowW(L"EDIT", L"", // Создаем поле ввода
                WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL, // Убираем ES_NUMBER, чтобы разрешить ввод любых чисел
                xOffset + j * 35, yOffset + i * 35, // Позиция
                30, 30, // Размеры
                hWnd, (HMENU)(100 + i * displayDimension + j), NULL, NULL); // Уникальный идентификатор
        }
    }

    // Проверяем, существует ли кнопка "Вычислить"
    if (!GetDlgItem(hWnd, 1)) { // Если кнопка не найдена
        // Создаем кнопку "Вычислить", если она не существует
        CreateWindowW(L"BUTTON", L"Вычислить", // Создаем кнопку
            WS_VISIBLE | WS_CHILD, // Стиль окна
            xOffset, yOffset + displayDimension * 30 + 60, // Позиция
            200, 40, // Размеры
            hWnd, (HMENU)(1), NULL, NULL); // Уникальный идентификатор
        // Создание кнопки "Вычислить из файла"
        CreateWindowW(L"BUTTON", L"Вычислить из файла", // Создаем кнопку
            WS_VISIBLE | WS_CHILD, // Стиль окна
            xOffset + 200, yOffset + displayDimension * 30 + 60, // Позиция
            200, 40, // Размеры
            hWnd, (HMENU)(2), NULL, NULL); // Уникальный идентификатор (2)
    }
}

// Функция для обновления полей ввода матрицы
void UpdateMatrixInputFields(HWND hWnd, SymQR* qr) {
    // Удаляем старые поля ввода
    for (int i = 0; i < 10; i++) { // Максимальная размерность 10
        for (int j = 0; j < 10; j++) {
            HWND hEdit = GetDlgItem(hWnd, 100 + i * 10 + j); // Получаем указатель на поле ввода
            if (hEdit) {
                DestroyWindow(hEdit); // Уничтожаем старые поля ввода
            }
        }
    }

    // Создаем новые поля ввода с новой размерностью
    CreateMatrixInputFields(hWnd, qr); // Создаем новые поля ввода
}

// Функция для отображения собственных значений
void ShowEigenvalues(HWND hWnd, SymQR* qr, HWND hResultTextBox) {
    string result = qr->displayEigenvalues(); // Получаем строку с собственными значениями
    // Отображение собственных значений в текстовом поле
    SetWindowTextA(hResultTextBox, result.c_str()); // Устанавливаем текст в текстовом поле
}

// Основная функция приложения
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR args, int ncmdshow) {
    // Регистрация класса окна
    WNDCLASS SoftwareMainClass = SymQRClass((HBRUSH)COLOR_WINDOW, LoadCursor(NULL, IDC_ARROW), hInst, LoadIcon(NULL, IDI_QUESTION),
        L"SymQR", SoftwareMainProcedure); // Создание класса окна

    if (!RegisterClassW(&SoftwareMainClass)) { return -1; } // Если регистрация не удалась, выходим
    MSG SoftwareMainMessage = { 0 }; // Инициализация структуры MSG

    // Создание главного окна
    CreateWindow(L"SymQR", L"QR Matrix", WS_OVERLAPPEDWINDOW | WS_VISIBLE, 10, 10, 800, 800, NULL, NULL, hInst, NULL);

    // Главный цикл обработки сообщений
    while (GetMessage(&SoftwareMainMessage, NULL, NULL, NULL)) { // Получаем сообщения
        TranslateMessage(&SoftwareMainMessage); // Переводим сообщения
        DispatchMessageW(&SoftwareMainMessage); // Обрабатываем сообщения
    }
    return 0; // Возвращаем 0 при завершении
}

// Функция обработки сообщений для главного окна
LRESULT CALLBACK SoftwareMainProcedure(HWND hWnd, UINT msg, WPARAM wp, LPARAM lp) {
    static SymQR* qr = nullptr; // Указатель на объект SymQR
    static HWND hResultTextBox = NULL; // Указатель на текстовое поле для вывода собственных значений
    static HWND hDimensionComboBox = NULL; // Указатель на выпадающий список для выбора размерности
    static HWND hDimensionLabel = NULL; // Указатель на статический текст для размерности
    static HWND hDimensionLabel1 = NULL; // Указатель на статический текст для матрицы
    int N = 100; // Начальная размерность матрицы

    switch (msg) { // Обработка сообщений
    case WM_CREATE: // Сообщение о создании окна
        qr = new SymQR(N); // Создаем объект SymQR для 10x10 матрицы
        CreateMatrixInputFields(hWnd, qr); // Создаем поля ввода для матрицы

        // Создание текстового поля для отображения собственных значений
        hResultTextBox = CreateWindowW(L"EDIT", L"",
            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_MULTILINE | ES_READONLY | WS_VSCROLL,
            20, 600, 440, 150, // Задаем координаты и размеры
            hWnd, NULL, NULL, NULL);

        // Создание текстового поля для отображения слова Матрица
        hDimensionLabel1 = CreateWindowW(L"STATIC", L"Матрица:",
            WS_CHILD | WS_VISIBLE,
            20, 10, 300, 20, // Задаем координаты и размеры
            hWnd, NULL, NULL, NULL);

        // Создание статического текста "Выберите размерность матрицы"
        hDimensionLabel = CreateWindowW(L"STATIC", L"Выберите размерность матрицы:",
            WS_CHILD | WS_VISIBLE,
            20, 470, 300, 20, // Задаем координаты и размеры
            hWnd, NULL, NULL, NULL);

        // Создание выпадающего списка для выбора размерности матрицы
        hDimensionComboBox = CreateWindowW(L"COMBOBOX", NULL,
            WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST,
            20, 500, 100, 100,
            hWnd, NULL, NULL, NULL);

        // Добавление элементов в выпадающий список
        for (int i = 2; i <= 10; i++) { // Добавляем размерности от 2 до 10
            wchar_t buffer[10]; // Буфер для хранения строки
            swprintf(buffer, 10, L"%d", i); // Форматируем число в строку
            SendMessage(hDimensionComboBox, CB_ADDSTRING, 0, (LPARAM)buffer); // Добавляем размерности
        }
        SendMessage(hDimensionComboBox, CB_SETCURSEL, 0, 0); // Установка начального значения
        break;

    case WM_COMMAND: // Обработка команд (например, нажатий кнопок)
        if (LOWORD(wp) == 1) { // ID кнопки "Вычислить"
            try {
                qr->readMatrix(hWnd); // Читаем матрицу из полей ввода
                qr->computeEigenvalues(); // Вычисляем собственные значения
                ShowEigenvalues(hWnd, qr, hResultTextBox); // Показываем собственные значения
            }
            catch (const runtime_error& e) { // Обработка исключений
                // Преобразование char* в wchar_t*
                wstring_convert<codecvt_utf8<wchar_t>, wchar_t> converter; // Конвертер для преобразования
                wstring wideError = converter.from_bytes(e.what()); // Преобразуем сообщение об ошибке

                MessageBox(hWnd, wideError.c_str(), L"Error", MB_OK | MB_ICONERROR); // Показываем сообщение об ошибке
            }
        }
        else if (LOWORD(wp) == 2) { // ID кнопки "Вычислить из файла"
            try {
                qr->readMatrixFromFile("file.txt"); // Читаем матрицу из файла
                qr->computeEigenvalues(); // Вычисляем собственные значения
                ShowEigenvalues(hWnd, qr, hResultTextBox); // Показываем собственные значения

                // Обновляем поля ввода с новой размерностью
                UpdateMatrixInputFields(hWnd, qr); // Обновляем поля ввода
            }
            catch (const runtime_error& e) {
                wstring_convert<codecvt_utf8<wchar_t>, wchar_t> converter;
                wstring wideError = converter.from_bytes(e.what());
                MessageBox(hWnd, wideError.c_str(), L"Error", MB_OK | MB_ICONERROR);
            }
        }
        else if (HIWORD(wp) == CBN_SELCHANGE && (HWND)lp == hDimensionComboBox) { // Если изменился выбор в комбобоксе
            int selectedDimension = SendMessage(hDimensionComboBox, CB_GETCURSEL, 0, 0) + 2; // Получаем выбранную размерность
            qr->resize(selectedDimension); // Изменяем размерность матрицы
            UpdateMatrixInputFields(hWnd, qr); // Обновляем поля ввода
        }
        break;

    case WM_DESTROY: // Сообщение о закрытии окна
        delete qr; // Освобождаем память
        PostQuitMessage(0); // Завершаем приложение
        break;
    default:
        return DefWindowProc(hWnd, msg, wp, lp); // Обработка стандартных сообщений
    }
    return 0; // Возвращаем 0
}