#include <Novice.h>
#include <cmath>
#include <limits>
#include <numbers>
#include <imgui.h>

const char kWindowTitle[] = "GC2A_12_ハラサワミツタカ";
const int kWindowWidth = 1280;
const int kWindowHeight = 720;;

using namespace std;

struct Vector3 {
	float x;
	float y;
	float z;
};
struct Matrix4x4 {
	float m[4][4];
};
struct Sphere {
	Vector3 center;
	float radius;
};
struct Segment {
	Vector3 origin;
	Vector3 diff;
};

//ヘッダーのincludeうまくいかないのでとりあえずmainに書く
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
}
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;
	return result;
}
Matrix4x4 Multiply(const Matrix4x4& a, const Matrix4x4& b)
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.m[i][j] = a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] + a.m[i][3] * b.m[3][j];
		}
	}
	return result;
}
float Determinant(const Matrix4x4 m)
{
	float det;
	det =
		m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2] -
		m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2] -

		m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2] +

		m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2] -
		m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2] -

		m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	return det;
}
Matrix4x4 Scalar(float a, const Matrix4x4& m)
{
	Matrix4x4 result;
	result.m[0][0] = a * (
		m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[1][3] * m.m[2][1] * m.m[3][2] -
		m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][3] -
		m.m[1][1] * m.m[2][3] * m.m[3][2]);
	result.m[0][1] = a * (
		-m.m[0][1] * m.m[2][2] * m.m[3][3] -
		m.m[0][2] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[2][1] * m.m[3][2] +
		m.m[0][3] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[2][1] * m.m[3][3] +
		m.m[0][1] * m.m[2][3] * m.m[3][2]);
	result.m[0][2] = a * (
		m.m[0][1] * m.m[1][2] * m.m[3][3] +
		m.m[0][2] * m.m[1][3] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[3][2] -
		m.m[0][3] * m.m[1][2] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][3] -
		m.m[0][1] * m.m[1][3] * m.m[3][2]);
	result.m[0][3] = a * (
		-m.m[0][1] * m.m[1][2] * m.m[2][3] -
		m.m[0][2] * m.m[1][3] * m.m[2][1] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] +
		m.m[0][3] * m.m[1][2] * m.m[2][1] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] +
		m.m[0][1] * m.m[1][3] * m.m[2][2]);

	result.m[1][0] = a * (
		-m.m[1][0] * m.m[2][2] * m.m[3][3] -
		m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[1][3] * m.m[2][0] * m.m[3][2] +
		m.m[1][3] * m.m[2][2] * m.m[3][0] +
		m.m[1][2] * m.m[2][0] * m.m[3][3] +
		m.m[1][0] * m.m[2][3] * m.m[3][2]);
	result.m[1][1] = a * (
		m.m[0][0] * m.m[2][2] * m.m[3][3] +
		m.m[0][2] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][2] -
		m.m[0][3] * m.m[2][2] * m.m[3][0] -
		m.m[0][2] * m.m[2][0] * m.m[3][3] -
		m.m[0][0] * m.m[2][3] * m.m[3][2]);
	result.m[1][2] = a * (
		-m.m[0][0] * m.m[1][2] * m.m[3][3] -
		m.m[0][2] * m.m[1][3] * m.m[3][0] -
		m.m[0][3] * m.m[1][0] * m.m[3][2] +
		m.m[0][3] * m.m[1][2] * m.m[3][0] +
		m.m[0][2] * m.m[1][0] * m.m[3][3] +
		m.m[0][0] * m.m[1][3] * m.m[3][2]);
	result.m[1][3] = a * (
		m.m[0][0] * m.m[1][2] * m.m[2][3] +
		m.m[0][2] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] -
		m.m[0][3] * m.m[1][2] * m.m[2][0] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] -
		m.m[0][0] * m.m[1][3] * m.m[2][2]);

	result.m[2][0] = a * (
		m.m[1][0] * m.m[2][1] * m.m[3][3] +
		m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][1] -
		m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][3] -
		m.m[1][0] * m.m[2][3] * m.m[3][1]);
	result.m[2][1] = a * (
		-m.m[0][0] * m.m[2][1] * m.m[3][3] -
		m.m[0][1] * m.m[2][3] * m.m[3][0] -
		m.m[0][3] * m.m[2][0] * m.m[3][1] +
		m.m[0][3] * m.m[2][1] * m.m[3][0] +
		m.m[0][1] * m.m[2][0] * m.m[3][3] +
		m.m[0][0] * m.m[2][3] * m.m[3][1]);
	result.m[2][2] = a * (
		m.m[0][0] * m.m[1][1] * m.m[3][3] +
		m.m[0][1] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][1] -
		m.m[0][3] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][3] -
		m.m[0][0] * m.m[1][3] * m.m[3][1]);
	result.m[2][3] = a * (
		-m.m[0][0] * m.m[1][1] * m.m[2][3] -
		m.m[0][1] * m.m[1][3] * m.m[2][0] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] +
		m.m[0][3] * m.m[1][1] * m.m[2][0] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] +
		m.m[0][0] * m.m[1][3] * m.m[2][1]);

	result.m[3][0] = a * (
		-m.m[1][0] * m.m[2][1] * m.m[3][2] -
		m.m[1][1] * m.m[2][2] * m.m[3][0] -
		m.m[1][2] * m.m[2][0] * m.m[3][1] +
		m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[1][1] * m.m[2][0] * m.m[3][2] +
		m.m[1][0] * m.m[2][2] * m.m[3][1]);
	result.m[3][1] = a * (
		m.m[0][0] * m.m[2][1] * m.m[3][2] +
		m.m[0][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][2] * m.m[2][0] * m.m[3][1] -
		m.m[0][2] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][2] -
		m.m[0][0] * m.m[2][2] * m.m[3][1]);
	result.m[3][2] = a * (
		-m.m[0][0] * m.m[1][1] * m.m[3][2] -
		m.m[0][1] * m.m[1][2] * m.m[3][0] -
		m.m[0][2] * m.m[1][0] * m.m[3][1] +
		m.m[0][2] * m.m[1][1] * m.m[3][0] +
		m.m[0][1] * m.m[1][0] * m.m[3][2] +
		m.m[0][0] * m.m[1][2] * m.m[3][1]);
	result.m[3][3] = a * (
		m.m[0][0] * m.m[1][1] * m.m[2][2] +
		m.m[0][1] * m.m[1][2] * m.m[2][0] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] -
		m.m[0][2] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][2] -
		m.m[0][0] * m.m[1][2] * m.m[2][1]);

	return result;
}
Matrix4x4 Inverse(const Matrix4x4& a)
{
	Matrix4x4 result;
	float det = 1.0f / Determinant(a);
	result = Scalar(det, a);

	return result;
}
Matrix4x4 Transpose(const Matrix4x4& a)
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.m[i][j] = a.m[j][i];
		}
	}
	return result;
}
Matrix4x4 MakeIdentity4x4()
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == j) {
				result.m[i][j] = 1.0f;
			}
			else {
				result.m[i][j] = 0;
			}
		}
	}
	return result;
}

Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result;

	result.x =
		vector.x * matrix.m[0][0] +
		vector.y * matrix.m[1][0] +
		vector.z * matrix.m[2][0] +
		1.0f * matrix.m[3][0];
	result.y =
		vector.x * matrix.m[0][1] +
		vector.y * matrix.m[1][1] +
		vector.z * matrix.m[2][1] +
		1.0f * matrix.m[3][1];
	result.z =
		vector.x * matrix.m[0][2] +
		vector.y * matrix.m[1][2] +
		vector.z * matrix.m[2][2] +
		1.0f * matrix.m[3][2];
	float w =
		vector.x * matrix.m[0][3] +
		vector.y * matrix.m[1][3] +
		vector.z * matrix.m[2][3] +
		1.0f * matrix.m[3][3];

	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

Matrix4x4 MakeScaleMatrix(const Vector3& scale)
{
	Matrix4x4 result;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i != j) {
				result.m[i][j] = 0;
			}
		}
	}
	result.m[0][0] = scale.x;
	result.m[1][1] = scale.y;
	result.m[2][2] = scale.z;
	result.m[3][3] = 1;
	return result;
}
Matrix4x4 MakeTranslateMatrix(const Vector3& translate)
{
	Matrix4x4 result;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == j) {
				result.m[i][j] = 1;
			}
			if (i != j) {
				result.m[i][j] = 0;
			}
		}
	}
	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	return result;
}
Matrix4x4 MakeRotateXMatrix(float radianX)
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.m[i][j] = 0;
		}
	}
	result.m[0][0] = 1;
	result.m[1][1] = std::cos(radianX);
	result.m[1][2] = std::sin(radianX);
	result.m[2][1] = -std::sin(radianX);
	result.m[2][2] = std::cos(radianX);
	result.m[3][3] = 1;

	return result;
}
Matrix4x4 MakeRotateYMatrix(float radianY)
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.m[i][j] = 0;
		}
	}
	result.m[0][0] = std::cos(radianY);
	result.m[0][2] = -std::sin(radianY);
	result.m[1][1] = 1;
	result.m[2][0] = std::sin(radianY);
	result.m[2][2] = std::cos(radianY);
	result.m[3][3] = 1;
	return result;
}
Matrix4x4 MakeRotateZMatrix(float radianZ)
{
	Matrix4x4 result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.m[i][j] = 0;
		}
	}
	result.m[0][0] = std::cos(radianZ);
	result.m[0][1] = std::sin(radianZ);
	result.m[1][0] = -std::sin(radianZ);
	result.m[1][1] = std::cos(radianZ);
	result.m[2][2] = 1;
	result.m[3][3] = 1;

	return result;
}
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate)
{
	Matrix4x4 result;

	Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix));
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);

	result = Multiply(Multiply(scaleMatrix, rotateXYZMatrix), translateMatrix);

	return result;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
	Matrix4x4 result = { 0 };
	if (nearClip == 0) {
		return result;
	}

	result.m[0][0] = (1.0f / aspectRatio) * 1.0f / std::tan(fovY / 2.0f);
	result.m[1][1] = 1.0f / std::tan(fovY / 2.0f);
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1.0f;
	result.m[3][2] = (-nearClip * farClip) / (farClip - nearClip);

	return result;
}
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
{
	Matrix4x4 result = { 0 };

	result.m[0][0] = 2.0f / (right - left);
	result.m[1][1] = 2.0f / (top - bottom);
	result.m[2][2] = 1.0f / (farClip - nearClip);
	result.m[3][0] = (left + right) / (left - right);
	result.m[3][1] = (top + bottom) / (bottom - top);
	result.m[3][2] = nearClip / (nearClip - farClip);
	result.m[3][3] = 1.0f;

	return result;
}
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
	Matrix4x4 result = { 0 };

	result.m[0][0] = width / 2.0f;
	result.m[1][1] = -(height / 2.0f);
	result.m[2][2] = maxDepth - minDepth;
	result.m[3][0] = left + width / 2.0f;
	result.m[3][1] = top + height / 2.0f;
	result.m[3][2] = minDepth;
	result.m[3][3] = 1.0f;

	return result;
}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 11;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

	Vector3 gridLineStart[kSubdivision + 1] = {};
	Vector3 gridLineEnd[kSubdivision + 1] = {};
	Vector3 ScreenStart[kSubdivision + 1] = {};
	Vector3 ScreenEnd[kSubdivision + 1] = {};

	//奥から手前に線を順に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		gridLineStart[xIndex].x = (-(float)kSubdivision / 2.0f + (float)xIndex) * kGridEvery;
		gridLineStart[xIndex].z = -kGridHalfWidth;
		gridLineEnd[xIndex].x = (-(float)kSubdivision / 2.0f + (float)xIndex) * kGridEvery;
		gridLineEnd[xIndex].z = kGridHalfWidth;

		ScreenStart[xIndex] = Transform(Transform(gridLineStart[xIndex], viewProjectionMatrix), viewportMatrix);
		ScreenEnd[xIndex] = Transform(Transform(gridLineEnd[xIndex], viewProjectionMatrix), viewportMatrix);

		Novice::DrawLine(int(ScreenStart[xIndex].x), int(ScreenStart[xIndex].y), int(ScreenEnd[xIndex].x), int(ScreenEnd[xIndex].y), WHITE);
	}
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		gridLineStart[zIndex].z = (-(float)kSubdivision / 2.0f + (float)zIndex) * kGridEvery;
		gridLineStart[zIndex].x = -kGridHalfWidth;
		gridLineEnd[zIndex].z = (-(float)kSubdivision / 2.0f + (float)zIndex) * kGridEvery;
		gridLineEnd[zIndex].x = kGridHalfWidth;

		ScreenStart[zIndex] = Transform(Transform(gridLineStart[zIndex], viewProjectionMatrix), viewportMatrix);
		ScreenEnd[zIndex] = Transform(Transform(gridLineEnd[zIndex], viewProjectionMatrix), viewportMatrix);

		Novice::DrawLine(int(ScreenStart[zIndex].x), int(ScreenStart[zIndex].y), int(ScreenEnd[zIndex].x), int(ScreenEnd[zIndex].y), WHITE);
	}
}
void DrawSphere(const Sphere sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	const uint32_t kSubdivision = 12;
	const float kLonEvery = 2.0f * numbers::pi_v<float> / kSubdivision;
	const float kLatEvery = numbers::pi_v<float> / kSubdivision;

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {

		float lat = numbers::pi_v<float> / 2.0f + kLatEvery * latIndex;

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {

			float lon = lonIndex * kLonEvery;
			Vector3 a, b, c;

			a = {
				{sphere.center.x + sphere.radius * cos(lat) * cos(lon)},
				{sphere.center.y + sphere.radius * sin(lat)},
				{sphere.center.z + sphere.radius * cos(lat) * sin(lon)}
			};
			b = {
				{sphere.center.x + sphere.radius * cos(lat + kLatEvery) * cos(lon)},
				{sphere.center.y + sphere.radius * sin(lat + kLatEvery)},
				{sphere.center.z + sphere.radius * cos(lat + kLatEvery) * sin(lon)}
			};
			c = {
				{sphere.center.x + sphere.radius * cos(lat) * cos(lon + kLonEvery)},
				{sphere.center.y + sphere.radius * sin(lat)},
				{sphere.center.z + sphere.radius * cos(lat) * sin(lon + kLonEvery)}
			};

			Vector3 screenA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);

			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenB.x), int(screenB.y), int(screenC.x), int(screenC.y), color);
		}
	}
}

Vector3 Project(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	float dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	float length = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;

	result.x = (dot / length) * v2.x;
	result.y = (dot / length) * v2.y;
	result.z = (dot / length) * v2.z;

	return result;
}
Vector3 ClosestPoint(const Vector3& point, const Segment segment) {
	Vector3 result;
	//ココが多分違う気がしないでもない
	result = Add(segment.origin, Project(Subtract(point, segment.origin), segment.diff));
	return result;
}



// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	//カメラ
	Vector3 cameraTranslate = { 0.0f,2.0f,-20.4f };
	Vector3 cameraRotate = { 0.0f,0.0f,0.0f };
	Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);

	Matrix4x4 viewMatrix = Inverse(cameraMatrix);
	Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
	Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

	//グリッド
	Vector3 gridScale = { 2.0f,2.0f,2.0f };
	Vector3 gridRotate = { 0.0f,0.0f,0.0f };
	Vector3 gridTranslate = { 0.0f,0.0f,0.0f };
	Matrix4x4 gridWorldMatrix = MakeAffineMatrix(gridScale, gridRotate, gridTranslate);
	Matrix4x4 gridWorldViewProjectionMatrix = Multiply(gridWorldMatrix, Multiply(viewMatrix, projectionMatrix));

	//球
	Vector3 sphereScale = { 1.0f,1.0f,1.0f };
	Vector3 sphereRotate = { 0.0f,0.0f,0.0f };
	Vector3 sphereTranslate = { 0.0f,1.0f,0.5f };
	Matrix4x4 sphereWorldMatrix = MakeAffineMatrix(sphereScale, sphereRotate, sphereTranslate);
	Matrix4x4 sphereWorldViewProjectionMatrix = Multiply(sphereWorldMatrix, Multiply(viewMatrix, projectionMatrix));
	Sphere sphere = { {0,0,0},0.5f };

	//02-00
	Segment segment = {
		{-2.0f,-1.0f,0.0f},
		{3.0f,2.0f,2.0f}
	};
	Vector3 point = { -1.5f,0.6f,0.6f };
	Vector3 project = Project(Subtract(point, segment.origin), segment.diff);
	Vector3 closestPoint = ClosestPoint(point, segment);

	Sphere pointSphere{ point,0.01f };
	Matrix4x4 pointWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, point);
	Matrix4x4 pointWorldViewProjectionMatrix = Multiply(pointWorldMatrix, Multiply(viewMatrix, projectionMatrix));

	Sphere closestPointSphere{ closestPoint,0.01f };
	Matrix4x4 closestPointWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, closestPoint);
	Matrix4x4 closestPointWorldViewProjectionMatrix = Multiply(closestPointWorldMatrix, Multiply(viewMatrix, projectionMatrix));

	Matrix4x4 segmentOriginWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, segment.origin);
	Matrix4x4 segmentOriginWorldViewProjectionMatrix = Multiply(segmentOriginWorldMatrix, Multiply(viewMatrix, projectionMatrix));
	Matrix4x4 segmentDiffWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, segment.diff);
	Matrix4x4 segmentDiffWorldViewProjectionMatrix = Multiply(segmentDiffWorldMatrix, Multiply(viewMatrix, projectionMatrix));
	Vector3 start = Transform(Transform(segment.origin, segmentOriginWorldViewProjectionMatrix), viewportMatrix);
	Vector3 end = Transform(Transform(Add(segment.origin,segment.diff), segmentDiffWorldViewProjectionMatrix),viewportMatrix);



	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///


		cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		viewMatrix = Inverse(cameraMatrix);
		projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);


		gridWorldMatrix = MakeAffineMatrix(gridScale, gridRotate, gridTranslate);
		gridWorldViewProjectionMatrix = Multiply(gridWorldMatrix, Multiply(viewMatrix, projectionMatrix));

		sphereWorldMatrix = MakeAffineMatrix(sphereScale, sphereRotate, sphereTranslate);
		sphereWorldViewProjectionMatrix = Multiply(sphereWorldMatrix, Multiply(viewMatrix, projectionMatrix));

		pointWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, point);
		pointWorldViewProjectionMatrix = Multiply(pointWorldMatrix, Multiply(viewMatrix, projectionMatrix));

		closestPointWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, closestPoint);
		closestPointWorldViewProjectionMatrix = Multiply(closestPointWorldMatrix, Multiply(viewMatrix, projectionMatrix));

		segmentOriginWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, segment.origin);
		segmentOriginWorldViewProjectionMatrix = Multiply(segmentOriginWorldMatrix, Multiply(viewMatrix, projectionMatrix));
		segmentDiffWorldMatrix = MakeAffineMatrix({ 1,1,1 }, { 0,0,0 }, segment.diff);
		segmentDiffWorldViewProjectionMatrix = Multiply(segmentDiffWorldMatrix, Multiply(viewMatrix, projectionMatrix));
		start = Transform(Transform(segment.origin, segmentOriginWorldViewProjectionMatrix), viewportMatrix);
		end = Transform(Transform(Add(segment.origin, segment.diff), segmentDiffWorldViewProjectionMatrix), viewportMatrix);



		ImGui::Begin("Window");
		ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
		ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
		ImGui::DragFloat3("SphereCenter", &sphere.center.x, 0.01f);
		ImGui::DragFloat3("SphereRotate", &sphereRotate.x, 0.01f);
		ImGui::DragFloat3("gridTranslate", &gridTranslate.x, 0.01f);
		ImGui::DragFloat("SphereRadius", &sphere.radius, 0.01f);
		ImGui::DragFloat3("project", &project.x,0.01f);

		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(gridWorldViewProjectionMatrix, viewportMatrix);
		DrawSphere(sphere, sphereWorldViewProjectionMatrix, viewportMatrix, WHITE);

		DrawSphere(pointSphere, pointWorldViewProjectionMatrix, viewportMatrix, RED);
		DrawSphere(closestPointSphere, closestPointWorldViewProjectionMatrix, viewportMatrix, GREEN);
		Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
