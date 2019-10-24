#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

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

class VecPoly
{
private:
	static vector<bool> my_xor(const vector<bool>& a, const vector<bool>& b, size_t offset = 1)
	{
		vector<bool> res;

		for (size_t i = offset; i < min(a.size(), b.size()); i++)
			res.push_back(a[i] != b[i]);

		return res;
	}

public:

	vector<bool> vec;

	VecPoly(size_t len) : vec(len, false) { }
	VecPoly(const vector<bool>& _vec) : vec(_vec) { }

	operator bool() const
	{
		for (bool b : vec)
			if (b)
				return true;

		return false;
	}

	VecPoly operator>>(int offset) const
	{
		vector<bool> res(vec.size(), 0);
		if (vec.size() >= offset)
		{
			copy(vec.begin() + offset, vec.end(), res.begin());
		}
		return res;
	}

	VecPoly operator<<(int offset) const
	{
		vector<bool> res(vec.size(), 0);
		if (vec.size() >= offset)
		{
			copy(vec.begin(), vec.end() - offset, res.begin() + offset);
		}
		return res;
	}

	VecPoly operator&(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < min(vec.size(), other.vec.size()); i++)
			res.push_back(vec[i] && other.vec[i]);

		return res;
	}

	VecPoly operator|(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < min(vec.size(), other.vec.size()); i++)
			res.push_back(vec[i] || other.vec[i]);

		return res;
	}

	VecPoly operator^(const VecPoly& other) const
	{
		return my_xor(vec, other.vec, 0);
	}

	VecPoly operator*(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < other.vec.size(); ++i)
		{
			for (size_t j = 0; j < vec.size(); ++j)
			{
				if (other.vec[i] && vec[j])
				{
					size_t offset = i + j;
					if (res.size() < offset + 1)
						res.resize(offset + 1, false);
					res[offset] = true;
				}
			}
		}

		return res;
	}

	VecPoly operator%(const VecPoly& divisor) const
	{
		size_t pick = divisor.vec.size();

		VecPoly tmp(pick);
		copy(vec.begin(), vec.begin() + min(pick, vec.size()), tmp.vec.begin());

		//cout << "begin: " << tmp.to_string() << endl;

		for (; pick < vec.size(); pick++)
		{
			//cout << tmp.to_string() << " " << pick << endl;
			if (tmp.vec[0])
			{
				tmp.vec = my_xor(divisor.vec, tmp.vec);
			}
			else
			{
				VecPoly t(pick);
				tmp.vec = my_xor(t.vec, tmp.vec);
			}
			tmp.vec.push_back(vec[pick]);
			//cout << tmp.to_string() << " " << pick+1 << endl;
		}

		if (tmp.vec[0])
		{
			tmp.vec = my_xor(divisor.vec, tmp.vec);
		}
		else
		{
			VecPoly t(pick);
			tmp.vec = my_xor(t.vec, tmp.vec);
		}

		return tmp;
	}

	string to_string() const
	{
		string res;

		for (bool b : vec)
			res += b ? "1" : "0";

		return res;
	}

	void crshift(int len)
	{
		*this = (*this >> len) | (*this << (int)(vec.size() - len));
	}

	void clshift(int len)
	{
		*this = (*this << len) | (*this >> (int)(vec.size() - len));
	}
};

class FireEncDec
{
private:
	// неприводимые полиномы g(x)
	static const map<int, vector<vector<bool>>> gxes;

	int n; // блоковая длина кода Файра
	int r; // длина проверочной группы в кодовой комбинации
	int payload_len; // длина полезной информации
	size_t gxf_msb; // старшая степень полинома g(x)f

	// полиномы для декодирования
	VecPoly gx;
	VecPoly xc;

	// полином для кодирования
	VecPoly gxf;

	static long filesize(const string& filename)
	{
		struct stat stat_buf;
		int rc = stat(filename.c_str(), &stat_buf);
		return rc == 0 ? stat_buf.st_size : -1;
	}

	void encode_block(VecPoly& info) const
	{
		info.vec.resize(info.vec.size() + gxf_msb, 0); // Производим увеличение степени информационного полинома I(X) на старшую степень образующего полинома G(X)Ф
		VecPoly remainder = info % gxf; // Выполняем деление расширенного информационного полинома Q(X) на образующий полином G(X)Ф
		info.vec.erase(info.vec.begin() + payload_len, info.vec.end());
		info.vec.insert(info.vec.end(), remainder.vec.begin(), remainder.vec.end());
	}

	void decode_block(VecPoly& cx) const
	{
		// Декодирование осуществляется путем деления кодового полинома С(Х) последовательно на обе компоненты образующего полинома кода Файра G(X)Ф: на (Х2+Х+1) и на (Х6+1).
		VecPoly info_remainder1 = cx % gx;
		VecPoly info_remainder2 = cx % xc;

		if (info_remainder1 && info_remainder2)
		{
			// обнаружена ошибка
			cout << "Found error!" << endl;

			VecPoly cx_orig = cx;
			size_t shift_count = 0;
			while (info_remainder1.vec != info_remainder2.vec && shift_count < n)
			{
				cx.crshift(1);
				info_remainder1 = cx % gx;
				info_remainder2 = cx % xc;
				shift_count++;

				// удалить нули в начале вектора
				while (!info_remainder1.vec.empty() && !info_remainder1.vec[0])
					info_remainder1.vec.erase(info_remainder1.vec.begin());

				while (!info_remainder2.vec.empty() && !info_remainder2.vec[0])
					info_remainder2.vec.erase(info_remainder2.vec.begin());
			}

			info_remainder1.vec.insert(info_remainder1.vec.begin(), cx_orig.vec.size() - info_remainder1.vec.size(), 0);
			info_remainder1.clshift(shift_count);
			cx = cx_orig ^ info_remainder1;
		}

		cx.vec.erase(cx.vec.end() - r, cx.vec.end()); // отбросить проверочный код и поделить на x^6 (на самом деле x^4)
	}

public:
	// P - длина исправляемого пакета
	FireEncDec(int P) : n(), r(), payload_len(), gxf_msb(), gx(0), xc(0), gxf(0)
	{
		int t = P; // старшая степень t не меньше длины исправляемого пакета P
		int c = 2 * P; // Старшая степень с второго полинома должна равняться или превышать удвоенную длину предполагаемого пакета ошибок Р
		int e = pow(2, t) - 1; // При этом с не должно делиться нацело на число е, где е = 2^t - 1
		while (c % e == 0)
			c++;

		gx = gxes.at(t)[0];

		xc.vec.resize(c + 1, false);
		xc.vec[0] = 1;
		xc.vec[c] = 1;

		gxf = gx * xc; // образующий полином кода Файра G(X)Ф определяется произведением двух примитивных полиномов: G(X)Ф = G(X)(Xс+1)
		gxf_msb = gxf.vec.size() - 1; // старшая степень полинома g(x)f

		n = e * c; // блоковая длина кода Файра будет являться произведением чисел e и с: n = е*с
		r = t + c; // длина проверочной группы r в кодовой комбинации будет являться суммой старших степеней обоих образующих полиномов: r = t+c
		payload_len = n - r; // длина полезной информации

		cout << "P = " << P << endl;
		cout << "t = " << t << endl;
		cout << "c = " << c << endl;
		cout << "e = " << e << endl;
		cout << "gx = " << gx.to_string() << endl;
		cout << "xc = " << xc.to_string() << endl;
		cout << "gxf = " << gxf.to_string() << endl;
		cout << "gxf_msb = " << gxf_msb << endl;
		cout << "n = " << n << endl;
		cout << "r = " << r << endl;
		cout << "n-r = " << payload_len << endl;
	}

	void encode(const string& filename, const string& out_filename)
	{
		ifstream file(filename, ios::binary);
		size_t total_bytes = filesize(filename);

		ofstream output_file(out_filename, ios::binary);
		output_file.write((char*)&total_bytes, sizeof(total_bytes)); // fixme endianess

		size_t in_block_size = payload_len; // in bits
		size_t out_block_size = n; // in bits
		size_t in_chunk_size = NOK(16, NOK(in_block_size, out_block_size)); // in bits
		size_t blocks_per_chunk = in_chunk_size / in_block_size;
		size_t out_chunk_size = out_block_size * blocks_per_chunk; // in bits
		size_t in_total_chunks = ceil(total_bytes * 8.0 / in_chunk_size);

#ifndef _DEBUG
#pragma omp parallel for
#endif
		for (int64_t chunk_num = 0; chunk_num < in_total_chunks; chunk_num++)
		{
			size_t in_byte_offset = chunk_num * (in_chunk_size / 8);
			vector<uint8_t> in_buf(in_chunk_size / 8);

			#pragma omp critical(FileIO)
			{
				file.clear();
				file.seekg(in_byte_offset);
				file.read((char*)in_buf.data(), in_buf.size());
			}

			vector<bool> in_chunk;
			in_chunk.reserve(in_chunk_size);
			for (uint8_t byte_data : in_buf)
			{
				for (int j = 7; j >= 0; j--)
					in_chunk.push_back((byte_data >> j) & 1);
			}

			vector<uint8_t> out_buf;
			int bits_wrote = 0;
			uint8_t out_byte = 0;
			int blocks = in_chunk.size() / in_block_size;
			for (int i = 0; i < blocks; i++)
			{
				VecPoly info(in_block_size);
				size_t in_bit_offset = i * in_block_size;
				copy(in_chunk.begin() + in_bit_offset, in_chunk.begin() + in_bit_offset + in_block_size, info.vec.begin());
				encode_block(info);

				for (bool b : info.vec)
				{
					out_byte |= b ? 1 : 0;
					bits_wrote++;

					if (bits_wrote == 8)
					{
						out_buf.push_back(out_byte);
						out_byte = 0;
						bits_wrote = 0;
					}

					out_byte <<= 1;
				}
			}

			size_t out_byte_offset = chunk_num * (out_chunk_size / 8) + sizeof(total_bytes);
			#pragma omp critical(FileIO)
			{
				output_file.seekp(out_byte_offset);
				output_file.write((char*) out_buf.data(), out_buf.size());
			}
		}
	}

	void decode(const string& filename, const string& out_filename)
	{
		ifstream file(filename, ios::binary);
		size_t out_total_bytes, total_bytes = filesize(filename) - sizeof(out_total_bytes);
		file.read((char*)&out_total_bytes, sizeof(out_total_bytes)); // fixme endianess

		ofstream output_file(out_filename, ios::binary);

		size_t in_block_size = n; // in bits
		size_t out_block_size = payload_len; // in bits
		size_t in_chunk_size = NOK(16, NOK(in_block_size, out_block_size)); // in bits
		size_t blocks_per_chunk = in_chunk_size / in_block_size;
		size_t out_chunk_size = out_block_size * blocks_per_chunk; // in bits
		size_t in_total_chunks = ceil(total_bytes * 8.0 / in_chunk_size);

#ifndef _DEBUG
#pragma omp parallel for
#endif
		for (int64_t chunk_num = 0; chunk_num < in_total_chunks; chunk_num++)
		{
			size_t in_byte_offset = chunk_num * (in_chunk_size / 8) + sizeof(out_total_bytes);
			vector<uint8_t> in_buf(in_chunk_size / 8);

			#pragma omp critical(FileIO)
			{
				file.clear();
				file.seekg(in_byte_offset);
				file.read((char*) in_buf.data(), in_buf.size());
			}

			vector<bool> in_chunk;
			in_chunk.reserve(in_chunk_size);
			for (uint8_t byte_data: in_buf)
			{
				for (int j = 7; j >= 0; j--)
					in_chunk.push_back((byte_data >> j) & 1);
			}

			size_t out_byte_offset = chunk_num * (out_chunk_size / 8);
			vector<uint8_t> out_buf;
			int bits_wrote = 0;
			uint8_t out_byte = 0;
			int blocks = in_chunk.size() / in_block_size;
			for (int i = 0; i < blocks; i++)
			{
				VecPoly info(in_block_size);
				size_t in_bit_offset = i * in_block_size;
				copy(in_chunk.begin() + in_bit_offset, in_chunk.begin() + in_bit_offset + in_block_size, info.vec.begin());
				decode_block(info);

				for (bool b : info.vec)
				{
					out_byte |= b ? 1 : 0;
					bits_wrote++;

					if (bits_wrote == 8 && out_byte_offset < out_total_bytes)
					{
						out_buf.push_back(out_byte);
						out_byte = 0;
						bits_wrote = 0;
						out_byte_offset++;
					}

					out_byte <<= 1;
				}
			}

			out_byte_offset = chunk_num * (out_chunk_size / 8);
			#pragma omp critical(FileIO)
			{
				output_file.seekp(out_byte_offset);
				output_file.write((char*) out_buf.data(), out_buf.size());
			}
		}
	}
};

const map<int, vector<vector<bool>>> FireEncDec::gxes = {
	{2, {{1, 1, 1}}},
	{3, {{1, 0, 1, 1}}},
	{4, {{1, 0, 0, 1, 1}, {1, 1, 1, 1, 1}}},
	{5, {{1, 0, 0, 1, 0, 1}}},
	{6, {{1, 0, 0, 0, 0, 1, 1}, {1, 0, 0, 1, 0, 0, 1}, {1, 1, 0, 1, 1, 1, 1}}},
	{7, {{1, 1, 1, 0, 0, 1, 1, 1}, {1, 0, 0, 0, 0, 0, 1, 1}, {1, 0, 0, 1, 1, 1, 0, 1}}}
};

int main(int argc, char** argv)
{
	string filename;
	int P;

	cout << "Filename: " << flush;
	cin >> filename;

	cout << "Expected error rate: " << flush;
	cin >> P;

	if (P < 2 || P > 7)
	{
		cout << "EER can be from 2 to 7" << endl;
		return 1;
	}

	FireEncDec fire(P);

	char mode;
	cout << "Mode (a - all, d - decode, e - encode): " << flush;
	cin >> mode;

	if (mode == 'a' || mode == 'e')
	{
		cout << "Encoding..." << flush;
		fire.encode(filename, filename + ".enc");
		filename += ".enc";
		cout << " done!" << endl;
	}

	if (mode == 'a' || mode == 'd')
	{
		cout << "Decoding..." << flush;
		fire.decode(filename, filename + ".dec");
		cout << " done!" << endl;
	}

	system("pause");
}
