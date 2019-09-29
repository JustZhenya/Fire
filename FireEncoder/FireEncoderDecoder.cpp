#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <bitset>

using namespace std;

// палиндром
typedef bitset<32> Poly;

//const size_t buf_size = 100 * 1024 * 1024;

// реверс битов в битсете
template<size_t N> void reverse(bitset<N>& b) {
	for (size_t i = 0; i < N / 2; ++i) {
		bool t = b[i];
		b[i] = b[N - i - 1];
		b[N - i - 1] = t;
	}
}

// Наибольший общий делитель
int NOD(int n1, int n2)
{
	int div;
	if (n1 == n2)  return n1;
	int d = n1 - n2;
	if (d < 0) {
		d = -d;
		div = NOD(n1, d);
	}
	else
		div = NOD(n2, d);
	return div;
}
// Наименьшее общее кратное
int NOK(int n1, int n2)
{
	return n1 * n2 / NOD(n1, n2);
}

// простое ли число
bool isPrime(int n)
{
	// Corner case 
	if (n <= 1)
		return false;

	// Check from 2 to n-1 
	for (int i = 2; i < n; i++)
		if (n % i == 0)
			return false;

	return true;
}

// палиндром с операторами сложения, деления, умножения, взятия остатка
class MyPoly
{
private:
	static void divide(const MyPoly& pa, const MyPoly& pb, Poly& presult, Poly& premainder)
	{
		uint32_t a_size = pa.get_msbit() + 1, b_size = pb.get_msbit() + 1;

		vector<bool> a, b, result;
		for (size_t i = 0; i < a_size; ++i)
			a.push_back(pa.p[i]);
		for (size_t i = 0; i < b_size; ++i)
			b.push_back(pb.p[i]);

		reverse(a.begin(), a.end());
		reverse(b.begin(), b.end());

		size_t append_zeroes = b.size() - 1;
		for (size_t i = 0; i < append_zeroes; ++i)
			a.push_back(false);

		while (b.size() <= a.size() && !a.empty())
		{
			if (a[0])
			{
				a.erase(a.begin());
				for (size_t j = 0; j < b.size() - 1; ++j)
					a[j] = (a[j] != b[j + 1]); // xor

				if (!a.empty())
					result.push_back(true);
			}
			else
			{
				a.erase(a.begin());
				result.push_back(false);
			}
		}

		reverse(result.begin(), result.end());
		reverse(a.begin(), a.end());

		presult = 0;
		for (size_t i = append_zeroes; i < result.size(); ++i)
			presult[i - append_zeroes] = result[i];

		premainder = 0;
		for (size_t i = 0; i < a.size(); ++i)
			premainder[i] = a[i];
	}

public:
	Poly p;

	MyPoly() = default;

	MyPoly(const Poly& _p): p(_p) { }

	MyPoly operator*(const MyPoly& other) const
	{
		MyPoly result;
		for (size_t i = 0; i < other.p.size(); ++i)
		{
			for (size_t j = 0; j < p.size(); ++j)
			{
				if(other.p[i] && p[j])
					result.p[i + j] = 1;
			}
		}
		return result;
	}

	// не уверен, правильно или нет
	MyPoly operator+(const MyPoly& other) const
	{
		return other.p | p;
	}

	MyPoly operator/(const MyPoly& divisor) const
	{
		Poly result, remainder;
		divide(p, divisor, result, remainder);
		return result;
	}

	MyPoly operator%(const MyPoly& divisor) const
	{
		Poly result, remainder;
		divide(p, divisor, result, remainder);
		return remainder;
	}

	/*MyPoly operator%(const MyPoly& divisor)
	{
		size_t a_size = get_msbit(), b_size = divisor.get_msbit();

		size_t result_idx = 0;
		Poly a = p, b = divisor.p, result;

		size_t append_zeroes = b_size - 1;
		a <<= append_zeroes;
		a_size += append_zeroes;

		reverse(a);
		reverse(b);

		while (b_size <= a_size && a_size)
		{
			if (a[a_size - 1])
			{
				//a >>= 1;
				--a_size;

				for (size_t j = 0; j > b_size - 1; --j)
					a[j] = (a[j] != b[j + 1]); // xor

				if (a_size)
					result[result_idx++] = true;
			}
			else
			{
				//a >>= 1;
				--a_size;

				if (a_size)
					result[result_idx++] = false;
			}
		}

		return result;
	}*/

	string to_string() const
	{
		string res;

		size_t msb = get_msbit() + 1;
		for (size_t i = 0; i < msb; ++i)
		{
			res = (p[i] ? "1" : "0") + res;
		}

		return res;
	}

	// возвращает старший бит
	// для 00110011 вернёт 5
	size_t get_msbit() const
	{
		size_t idx = 0;
		for (size_t i = 0; i < p.size(); ++i)
		{
			if (p[i])
				idx = i;
		}
		return idx;
	}
};

class FireEncDec
{
private:

	// неприводимые полиномы g(x)
	static const map<int, vector<Poly>> gxes;

	int P; // длина исправляемого пакета
	int t; // старшая степень полинома gx
	int c; // Старшая степень второго полинома
	int e; // При этом с не должно делиться нацело на число е, где е = 2^t - 1
	int n; // При этом блоковая длина кода Файра будет являться произведением чисел e и с: n = е*с
	int r; // а длина проверочной группы r в кодовой комбинации будет являться суммой старших степеней обоих образующих полиномов: r = t+c
	int payload_len; // длина полезной информации

	MyPoly gx;
	MyPoly xc;
	MyPoly gxf;
	MyPoly msb_poly;

public:
	FireEncDec(int _P): P(_P), t(0), c(0), e(0), n(0), r(0), payload_len(0), gx(0), xc(0), gxf(0), msb_poly(0)
	{
		// старшая степень t не меньше длины исправляемого пакета P
		t = P;
		gx = gxes.at(t)[0];

		// Старшая степень с второго полинома должна равняться или превышать удвоенную длину предполагаемого пакета ошибок Р:
		c = 2 * P;
		// При этом с не должно делиться нацело на число е, где е = 2^t - 1
		e = pow(2, t) - 1;
		while (c % e == 0)
			c++;

		xc.p[0] = 1;
		xc.p[c] = 1;

		// При этом блоковая длина кода Файра будет являться произведением чисел e и с: n = е*с
		n = e * c; // длина кода Файра

		// а длина проверочной группы r в кодовой комбинации будет являться суммой старших степеней обоих образующих полиномов: r = t+c
		r = t + c; // длина избыточной информации

		payload_len = n - r; // длина полезной информации

		// образующий полином кода Файра G(X)Ф определяется произведением двух примитивных полиномов: G(X)Ф = G(X)(Xс+1)
		gxf = gx * xc;

		// ...старшую степень образующего полинома G(X)Ф
		msb_poly.p[gxf.get_msbit()] = 1;

		cout << "P = " << P << endl;
		cout << "t = " << t << endl;
		cout << "c = " << c << endl;
		cout << "e = " << e << endl;
		cout << "n = " << n << endl;
		cout << "r = " << r << endl;
		cout << "n-r = " << payload_len << endl;
		cout << "gx = " << gx.to_string() << endl;
		cout << "xc = " << xc.to_string() << endl;
		cout << "gxf = " << gxf.to_string() << endl;
		cout << "msb_poly = " << msb_poly.to_string() << endl;
	}

	void encode(const string& filename)
	{
		MyPoly info;

		ifstream file(filename, ios::binary);
		ofstream output_file(filename + ".enc", ios::binary);
		char in_byte = 0, out_byte = 0;
		int bits_read = 0, bits_wrote = 0;
		while (file.get(in_byte))
		{
			for (int i = 7; i >= 0; i--)
			{
				bool bit = (in_byte >> i) & 1;
				info.p |= bit;
				bits_read++;
				if (bits_read == payload_len)
				{
					// Производим увеличение степени информационного полинома I(X) на старшую степень образующего полинома G(X)Ф
					MyPoly qx = info * msb_poly;

					// Кодовый полином C(X) получается путем добавления остатка от деления Q(X) на G(X)Ф к расширенному информационному полиному Q(X).
					MyPoly cx = (qx % gxf) + qx;

					for (int j = n - 1; j >= 0; j--)
					{
						// write bit
						out_byte |= cx.p[j];
						bits_wrote++;

						if (bits_wrote == 8)
						{
							output_file.put(out_byte);
							out_byte = 0;
							bits_wrote = 0;
						}

						out_byte <<= 1;
					}

					bits_read = 0;
					info.p = 0;
				}
				info.p <<= 1;
			}
		}
	}

	void decode(const string& filename)
	{
		MyPoly cx;

		ifstream file(filename, ios::binary);
		ofstream output_file(filename + ".dec", ios::binary);
		char in_byte = 0, out_byte = 0;
		int bits_read = 0, bits_wrote = 0;
		while (file.get(in_byte))
		{
			for (int i = 7; i >= 0; i--)
			{
				bool bit = (in_byte >> i) & 1;
				cx.p |= bit;
				bits_read++;
				if (bits_read == n)
				{
					// Декодирование осуществляется путем деления кодового полинома С(Х) последовательно на обе компоненты образующего полинома кода Файра G(X)Ф: на (Х2+Х+1) и на (Х6+1).
					MyPoly info = cx / gx;
					info = info / xc;

					for (int j = payload_len - 1; j >= 0; j--)
					{
						out_byte |= info.p[j];
						bits_wrote++;

						if (bits_wrote == 8)
						{
							output_file.put(out_byte);
							out_byte = 0;
							bits_wrote = 0;
						}

						out_byte <<= 1;
					}

					bits_read = 0;
					cx.p = 0;
				}
				cx.p <<= 1;
			}
		}
	}
};

const map<int, vector<Poly>> FireEncDec::gxes = {
	{2, {0b111}},
	{3, {0b1011}},
	{4, {0b10011, 0b11111}},
	{5, {0b100101}},
	{6, {0b1000011, 0b1001001, 0b1101111}},
	{7, {0b11100111, 0b10000011, 0b10011101}}
};

int main()
{
	string filename;
	int P;

	cout << "Filename: ";
	cin >> filename;

	cout << "Expected error rate: ";
	cin >> P;

	if (P < 1 || P > 8)
	{
		cout << "EER can be from 1 to 8" << endl;
		return 1;
	}

	FireEncDec fire(P);

	cout << "Encoding..." << flush;
	fire.encode(filename);
	cout << " done!" << endl;

	cout << "Decoding..." << flush;
	fire.decode(filename + ".enc");
	cout << " done!" << endl;

	system("pause");

	return 0;
}
