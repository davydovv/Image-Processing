#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <windows.h>
#include <vector>
#include <fstream>

using namespace std;

struct rgb 
{
	unsigned char b;
	unsigned char g;
	unsigned char r;
};

struct ycbcr 
{
	unsigned char y;
	unsigned char cb;
	unsigned char cr;
};



const unsigned colorMaskRed = 0x1;
const unsigned colorMaskGreen = 0x2;
const unsigned colorMaskBlue = 0x4;
const unsigned colorMaskY = 0x8;
const unsigned colorMaskCb = 0x10;
const unsigned colorMaskCr = 0x20;

void only_red(rgb *pix) {
	pix->g = 0;
	pix->b = 0;
}

void only_blue(rgb *pix) {
	pix->g = 0;
	pix->r = 0;
}

void only_green(rgb *pix) {
	pix->r = 0;
	pix->b = 0;
}

class BmpFile
{
	bool             is_valid = false;
	BITMAPFILEHEADER file_hdr;
	BITMAPINFOHEADER file_info_hdr;
	vector<char>     file_other_hdrs;
	vector<rgb>      file_data;
	vector<char>     file_tail;

	rgb              dummy_pixel;

public:
	size_t getPixelCount();
	rgb&   getPixel(size_t i);
	size_t getW();
	size_t getH();
	void   load(string filename);
	bool   save(string filename);
};

class IPixelComponentBlock
{
public:
	virtual unsigned char get(size_t x, size_t y) = 0;
	virtual void set(size_t x, size_t y, char value) = 0;
	virtual size_t get_size_x() = 0;
	virtual size_t get_size_y() = 0;
	void print()
	{
		for (int k = 0; k < get_size_x(); k++) {
			for (int l = 0; l < get_size_y(); l++)
				printf("%03d ", get(k, l));
			printf("\n");
		}
		printf("\n");
	}
};

class BmpFilePixelComponentBlock : public IPixelComponentBlock
{
	BmpFile *data;
	size_t offset_x;
	size_t offset_y;
	size_t size_x;
	size_t size_y;
	int component_type;
public:
	BmpFilePixelComponentBlock(BmpFile *data, size_t offset_x, size_t offset_y, int component_type, int sizex = 0, int sizey = 0)
	{
		this->data = data;
		this->offset_x = offset_x;
		this->offset_y = offset_y;
		this->size_x = sizex == 0 ? data->getW() - offset_x : sizex;
		this->size_y = sizey == 0 ? data->getH() - offset_y : sizey;
		this->component_type = component_type;
	}

	unsigned char get(size_t x, size_t y)
	{
		rgb& pixel = data->getPixel((offset_y + y) * data->getW() + offset_x + x);
		if (component_type == colorMaskRed) return pixel.r;
		else if (component_type == colorMaskGreen) return pixel.g;
		else return pixel.b;
	}

	void set(size_t x, size_t y, char value)
	{
		rgb& pixel = data->getPixel((offset_y + y) * data->getW() + offset_x + x);
		pixel.r = pixel.g = pixel.b = value;
	}

	virtual size_t get_size_x() { return size_x; }
	virtual size_t get_size_y() { return size_y; }
};

class ArrayComponentBlock : public IPixelComponentBlock
{
	char *data;
	size_t w;
	size_t offset_x;
	size_t offset_y;
public:
	ArrayComponentBlock(char *data, size_t w, size_t offset_x, size_t offset_y)
	{
		this->data = data;
		this->w = w;
		this->offset_x = offset_x;
		this->offset_y = offset_y;
	}

	unsigned char get(size_t x, size_t y)
	{
		return data[(offset_y + y) * w + offset_x + x];
	}
	void set(size_t x, size_t y, double value)
	{
		data[(offset_y + y) * w + offset_x + x] = value;
	}
};



class BmpProcessor
{
	ycbcr rgb2ycbcr(rgb src);
	void extractColors(BmpFile &data, unsigned colorMask);

public:
	double PSNR(BmpFile &data_in, BmpFile &data_out, unsigned char color);
	
	void extractGreen(BmpFile &data) { extractColors(data, colorMaskGreen); }
	
	void extractBlue(BmpFile &data)  { extractColors(data, colorMaskBlue); }
	
	void extractRed(BmpFile &data)   { extractColors(data, colorMaskRed); }
	
	void extractYCbCr(BmpFile &data, unsigned color);
	
	unsigned char getYCbCr(BmpFile &data, unsigned type);
	
	void extractY(BmpFile &data)     { extractYCbCr(data, colorMaskY); }
	
	void extractCb(BmpFile &data)    { extractYCbCr(data, colorMaskCb); }
	
	void extractCr(BmpFile &data)    { extractYCbCr(data, colorMaskCr); }
	
	double MO(BmpFile &data, unsigned color);
	
	double MO(IPixelComponentBlock &data); // mean value
	
	double sigma(BmpFile &data, unsigned color);
	double sigma(IPixelComponentBlock &data);
	
	double BmpProcessor::Korr( BmpFile &data, unsigned color1, unsigned color2, int x, int y); //correlation
	
	void mergeYCbCr(BmpFile &dst,  BmpFile &ydata,  BmpFile &cbdata,  BmpFile &crdata);

	void decimateByDublication(BmpFile &data, int k);
	void decimateByAvg(BmpFile &data, int k);

	void mirror_vert(BmpFile &data);
	void mirror_goriz(BmpFile &data);
	
	void turn_90(BmpFile &data);
	void turn_180(BmpFile &data);
	void turn_270(BmpFile &data);

	vector<int> DPCM(IPixelComponentBlock &data, int r);

};

ycbcr BmpProcessor::rgb2ycbcr(rgb src)
{
	ycbcr res;
	res.y = 0.5 + 0.299 * src.r + 0.587 * src.g + 0.114 * src.b;
	res.cb = 128.5 + 0.5643 * (src.b - res.y);
	res.cr = 128.5 + 0.7132 * (src.r - res.y);
	return res;
}

void BmpProcessor::extractColors(BmpFile &data, unsigned colorMask)
{
	for (size_t i = 0; i < data.getPixelCount(); i++)
	{
		rgb &t = data.getPixel(i);
		if ((colorMask & colorMaskRed) == 0)
			t.r = 0;
		if ((colorMask & colorMaskGreen) == 0)
			t.g = 0;
		if ((colorMask & colorMaskBlue) == 0)
			t.b = 0;
	}
}

void BmpProcessor::extractYCbCr(BmpFile &data, unsigned type) 
{
	for (size_t i = 0; i < data.getPixelCount(); i++)
	{
		rgb &t = data.getPixel(i);
		ycbcr pix = rgb2ycbcr(t);
		if (colorMaskY == type) t.r = t.g = t.b = pix.y;
		else if (colorMaskCb == type) t.r = t.g = t.b = pix.cb;
		else t.r = t.g = t.b = pix.cr;
	}
}

unsigned char BmpProcessor::getYCbCr(BmpFile &data, unsigned type)
{
	for (size_t i = 0; i < data.getPixelCount(); i++)
	{
		rgb &t = data.getPixel(i);
		ycbcr pix = rgb2ycbcr(t);
		if (colorMaskY == type) return pix.y;
		else if (colorMaskCb == type) return pix.cb;
		else return pix.cr;
	}
}

void BmpProcessor::mergeYCbCr(BmpFile &dst, BmpFile &ydata,  BmpFile &cbdata,  BmpFile &crdata)
{
	for (size_t i = 0; i < ydata.getPixelCount(); i++)
	{
		unsigned char y = ydata.getPixel(i).r;
		unsigned char cb = cbdata.getPixel(i).r;
		unsigned char cr = crdata.getPixel(i).r;

		rgb &t = dst.getPixel(i);

		t.g = y - 0.714 * (cr - 128) - 0.334 * (cb - 128);
		t.r = y + 1.402	* (cr - 128);
		t.b = y + 1.772	* (cb - 128);


	}
}
double BmpProcessor::PSNR(BmpFile &data_in, BmpFile &data_out, unsigned char color)
{
	double k = 0;
	for (int i = 0; i < data_in.getH()*data_in.getW(); i++) {
		if (color == colorMaskRed)
			k += pow((int)data_in.getPixel(i).r - data_out.getPixel(i).r, 2);
		else if (color == colorMaskBlue)
			k += pow((int)data_in.getPixel(i).b - data_out.getPixel(i).b, 2);
		else
			k += pow((int)data_in.getPixel(i).g - data_out.getPixel(i).g, 2);
	}

	return 10 * log10(data_in.getH()*data_in.getW()*pow(255, 2) / k);
}

double BmpProcessor::Korr(BmpFile &data, unsigned color1, unsigned color2, int x, int y)
{
	BmpFilePixelComponentBlock block1(&data, 0, 0, color1, data.getW() - x, data.getH() - y);
	BmpFilePixelComponentBlock block2(&data, x, y, color2, data.getW() - x, data.getH() - y);
	double math1 = MO(block1);
	double math2 = MO(block2);

	double res = 0;
	for (int i = 0; i < data.getW() - x; i++) {
		for (int j = 0; j < data.getH() - y; j++)
			res += (block1.get(i, j) - math1) * (block2.get(i, j) - math2);
	}

	res /= block1.get_size_x() * block1.get_size_y();
	res /= sigma(block1) * sigma(block2);
	return res;
}

	
double BmpProcessor::MO(BmpFile &data, unsigned color)
{
	double res = 0;
	for (size_t i = 0; i < data.getPixelCount(); i++)
	{
		if (color == colorMaskRed) res += data.getPixel(i).r;
		else  if (color == colorMaskBlue) res += data.getPixel(i).b;
		else res += data.getPixel(i).g;
	}
	res /= data.getPixelCount();
	return res;
}

double BmpProcessor::MO(IPixelComponentBlock &data)
{
	double res = 0;
	for (int i = 0; i < data.get_size_x(); i++) {
		for (int j = 0; j < data.get_size_y(); j++) 
			res += data.get(i,j);
	}

	res /= data.get_size_x() * data.get_size_y();
	return res;
}

double BmpProcessor::sigma(BmpFile &data, unsigned color)
{
	double res = 0;
	double mo = MO(data, color);
	for (size_t i = 0; i < data.getPixelCount(); i++)
	{
		double t;
		if (color == colorMaskRed) t = data.getPixel(i).r - mo;
		else  if (color == colorMaskBlue) t = data.getPixel(i).b - mo;
		else t = data.getPixel(i).g - mo;
		res += t * t;
	}
	res /= data.getPixelCount() - 1;
	return sqrt(res);
}

double BmpProcessor::sigma(IPixelComponentBlock &data)
{
	double res = 0;
	double mo = MO(data);
	for (int i = 0; i < data.get_size_x(); i++) {
		for (int j = 0; j < data.get_size_y(); j++) {
			res += pow(data.get(i, j) - mo, 2);
		}
	}
	res /= (data.get_size_x() * data.get_size_y() - 1);
	return sqrt(res);
}

void BmpProcessor::decimateByDublication(BmpFile &data, int k)
{
	for (int i = 0; i < data.getW() - k; i += k) {
		for (int j = 0; j < data.getH() - k; j += k) {
			BmpFilePixelComponentBlock block(&data, i, j, 1);
			for (int l = 0; l <= k - 1; l++)
			for (int m = 0; m <= k - 1; m++)
				block.set(l, m, block.get(k - 1, k - 1));
		}
	}
}

void BmpProcessor::decimateByAvg(BmpFile &data, int k)
{
	for (int i = 0; i < data.getW() - k; i += k) {
		for (int j = 0; j < data.getH() - k; j += k) {
			BmpFilePixelComponentBlock block(&data, i, j, 1);
			double s = 0;
			for (int l = 0; l <= k - 1; l++)
			for (int m = 0; m <= k - 1; m++)
				s+=block.get(l, m);
			s /= k * k;
			for (int l = 0; l <= k - 1; l++)
			for (int m = 0; m <= k - 1; m++)
				block.set(l, m, s);
		}
	}
}

void BmpProcessor::mirror_vert(BmpFile &data)
{
	for (int i = 0; i < data.getW() / 2; i++) {
		for (int j = 0; j < data.getH(); j++) {
			rgb tmp = data.getPixel(j*data.getW() + i);
			data.getPixel(j * data.getW() + i) = data.getPixel(j*data.getW() + data.getW() - i - 1);
			data.getPixel(j * data.getW() + data.getW() - i - 1) = tmp;
		}
	}
}

void BmpProcessor::mirror_goriz(BmpFile &data)
{
	for (int i = 0; i < data.getW(); i++) {
		for (int j = 0; j < data.getH() / 2; j++) {
			data.getPixel(j*data.getW() + i) = data.getPixel((data.getH() - j - 1) * data.getW() + i);
		}
	}
}

void BmpProcessor::turn_90(BmpFile &data)
{
	for (int i = 0; i < data.getW(); i++) {
		for (int j = i; j < data.getH(); j++) {
			rgb tmp = data.getPixel(j*data.getW() + i);
			data.getPixel(j * data.getW() + i) = data.getPixel(i*data.getW() + j);
			data.getPixel(i*data.getW() + j) = tmp;
		}
	}
	mirror_vert(data);
}

vector<int> BmpProcessor::DPCM(IPixelComponentBlock &data, int r)
{
	vector<int> res;
	res.resize((data.get_size_x() - 1) * (data.get_size_y() - 1));
	for (int i = 1; i < data.get_size_y(); i++)
	{
		for (int j = 1; j < data.get_size_x(); j++)
		{
			int t = data.get(j, i);
			switch (r)
			{
			case 1: t -= data.get(j - 1, i); break;
			case 2: t -= data.get(j, i - 1); break;
			case 3: t -= data.get(j - 1, i - 1); break;
			case 4: t -= (data.get(j - 1, i - 1) + data.get(j - 1, i) + data.get(j, i - 1)) / 3; break;
			}
			res.at((i -1) * (data.get_size_x() - 1) + j - 1) = t;
		}
	}
	return res;
}

size_t BmpFile::getPixelCount()
{
	if (!is_valid)
		return 0;
	return file_info_hdr.biWidth * file_info_hdr.biHeight;
}

rgb& BmpFile::getPixel(size_t i)
{
	if (!is_valid || i >= file_info_hdr.biWidth * file_info_hdr.biHeight)
		return dummy_pixel;
	return file_data.at(i);
}


size_t BmpFile::getW()
{
	if (!is_valid)
		return 0;
	return file_info_hdr.biWidth;

}

size_t BmpFile::getH()
{
	if (!is_valid)
		return 0;
	return file_info_hdr.biHeight;

}

void BmpFile::load(string filename)
{
	is_valid = false;
	FILE *f;
	cout << "Loading file " << filename.data() << endl;
	fopen_s(&f, filename.data(), "rb");
	if (f == NULL)
		return;
	cout << "\tReading header" << endl;
	if (!fread(&file_hdr, sizeof(file_hdr), 1, f)) { fclose(f); return; }

	cout << "\tReading info header" << endl;
	if (!fread(&file_info_hdr, sizeof(file_info_hdr), 1, f)) { fclose(f); return; }

	int count_pix = file_info_hdr.biWidth * file_info_hdr.biHeight;

	cout << "\tReading pixels data" << endl;
	file_data.resize(count_pix);
	fseek(f, file_hdr.bfOffBits, SEEK_SET);
	if (fread(file_data.data(), sizeof(rgb), file_data.size(), f) != count_pix) { fclose(f); return; }

	file_other_hdrs.resize(file_info_hdr.biSize - sizeof(file_info_hdr));
	if (file_other_hdrs.size() > 0)
	{
		cout << "\tReading extra headers" << endl;
		fseek(f, sizeof(file_hdr)+sizeof(file_info_hdr), SEEK_SET);
		if (fread(file_other_hdrs.data(), 1, file_other_hdrs.size(), f) != file_other_hdrs.size()) { fclose(f); return; }
	}

	int src_dataend = file_info_hdr.biSize + sizeof(file_hdr)+count_pix * sizeof(rgb);
	file_tail.resize(file_hdr.bfSize - src_dataend);
	if (file_tail.size() > 0)
	{
		cout << "\tReading tail data" << endl;
		fseek(f, src_dataend, SEEK_SET);
		if (fread(file_tail.data(), 1, file_tail.size(), f) != file_tail.size()) { fclose(f); return; }
	}
	fclose(f);
	cout << "Loading finished. OK" << endl;
	is_valid = true;
}

bool BmpFile::save(string filename)
{
	if (!is_valid)
		return false;
	FILE * f;
	fopen_s(&f, filename.data(), "wb");

	fwrite(&file_hdr, sizeof(file_hdr), 1, f);
	fwrite(&file_info_hdr, sizeof(file_info_hdr), 1, f);
	fwrite(file_other_hdrs.data(), file_other_hdrs.size(), 1, f);
	fwrite(file_data.data(), file_data.size() * sizeof(rgb), 1, f);
	fwrite(file_tail.data(), file_tail.size(), 1, f);


	fclose(f);
	return true;
}


int main()
{
	//examples of usage
	BmpFile bf, bt, bcr, by, bcb, bx, btcb, btcr;
	BmpProcessor bp;
	bf.load("output.bmp");

	bt = bf; bp.extractGreen(bt); bt.save("lena_g.bmp");
	bt = bf; bp.extractRed(bt); bt.save("lena_r.bmp");
	bt = bf; bp.extractBlue(bt); bt.save("lena_b.bmp");

	bt = bf; bp.extractY(bt); bt.save("lena_y.bmp");
	bt = bf; bp.extractCb(bt); bt.save("lena_cb.bmp");
	bt = bf; bp.extractCr(bt); bt.save("lena_cr.bmp");
	by.load("lena_y.bmp");
	bcb.load("lena_cb.bmp");
	bcr.load("lena_cr.bmp");
	bp.mergeYCbCr(bt, by, bcb, bcr);

	bt.save("lena_new.bmp");

	bt = bf;
	bp.mirror_vert(bt);
	bt.save("lena_mirror.bmp");

// ...etc
	return 0;
	}

