#include "SettingData.h"
#include <fstream>
#include <numbers>

namespace png {
	LoadData::LoadData(std::string jsonName) {
		std::ifstream jsonStream(jsonName);
		nlohmann::json json;
		jsonStream >> json;

		from_json(json, data);
	}

	void LoadData::SaveSampleJson(std::string fileName) {
		SettingData tmp;
		tmp.width = 960;
		tmp.height = 540;
		tmp.samples = 10;
		tmp.superSamples = 4;
		tmp.camera.origin = vec3(0, 0, 0);
		tmp.camera.target = vec3(0, 0, 1);
		tmp.camera.upVec = vec3(0, 1, 0);
		tmp.camera.fov = 60;

		/*
		tmp.object.push_back({
			vec3(0,0,10)
			, 1.0
			, Material{
				vec3(1.0,0.2,0.2)
				,vec3(0,0,0)
			}
			});

		//light
		tmp.object.push_back({
			vec3(0,11,10)
			, 10.0
			, Material{
				vec3(1.0,1.0,1.0)
				,vec3(1,1,1)
			}
			});

		//wall
		{
			float size = 99999;
			float space = 10;
			tmp.object.push_back({
				vec3(size+0.5*space,0,0)
				, size
				, Material{
					vec3(1.0,0.5,0.5)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(-(size + 0.5 * space),0,0)
				, size
				, Material{
					vec3(0.5,1.0,0.5)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,size + 0.5 * space,0)
				, size
				, Material{
					vec3(0.5,0.5,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,-(size + 0.5 * space),0)
				, size
				, Material{
					vec3(0.5,1.0,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,size + 0.5 * space)
				, size
				, Material{
					vec3(1.0,0.5,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,-(size + 0.5 * space))
				, size
				, Material{
					vec3(1.0,1.0,0.5)
					,vec3()
				}
				});
		}
		*/

		nlohmann::json tmpJson;
		to_json(tmpJson, tmp);
		std::ofstream ostream(fileName);
		ostream << tmpJson;
	}

	void LoadData::to_json(nlohmann::json& json, const SettingData& data) {
		json = {
			{"00 width", data.width}
			,{"00 height", data.height}
			,{"00 samples",data.samples}
			,{"00 superSamples",data.superSamples}
			,{"01 camera",{
				{"origin",{data.camera.origin.x,data.camera.origin.y,data.camera.origin.z}}
				,{"target",{data.camera.target.x,data.camera.target.y,data.camera.target.z}}
				,{"upVec",{data.camera.upVec.x,data.camera.upVec.y,data.camera.upVec.z}}
				,{"fov",data.camera.fov}
			}}
		};
		for (int i = 0; i < data.object.size(); ++i) {
			auto& obj = data.object[i];
			/*
			json["02 scene"]["00 object"][i]["00 position"] = { obj.position.x,obj.position.y,obj.position.z };
			json["02 scene"]["00 object"][i]["01 size"] = obj.size;
			json["02 scene"]["00 object"][i]["02 material"]["color"] = { obj.material.color.x ,obj.material.color.y,obj.material.color.z };
			json["02 scene"]["00 object"][i]["02 material"]["emission"] = { obj.material.emission.x,obj.material.emission.y,obj.material.emission.z };
			*/
		}
	}
	void LoadData::from_json(const nlohmann::json& json, SettingData& data) {
		for (auto& it : json.items()) {
			if (it.key() == "00 width") data.width = it.value();
			else if (it.key() == "00 height") data.height = it.value();
			else if (it.key() == "00 samples") data.samples = it.value();
			else if (it.key() == "00 superSamples") data.superSamples = it.value();
			else if (it.key() == "01 camera") {
				for (auto& it_cam : it.value().items()) {
					if (it_cam.key() == "origin") {
						data.camera.origin.x = it_cam.value()[0];
						data.camera.origin.y = it_cam.value()[1];
						data.camera.origin.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "target") {
						data.camera.target.x = it_cam.value()[0];
						data.camera.target.y = it_cam.value()[1];
						data.camera.target.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "upVec") {
						data.camera.upVec.x = it_cam.value()[0];
						data.camera.upVec.y = it_cam.value()[1];
						data.camera.upVec.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "fov") {
						data.camera.fov = it_cam.value();
					}
				}
			}
			else if (it.key() == "02 scene") {
				for (auto& it_scene : it.value().items()) {
					if (it_scene.key() == "00 object") {
						for (auto& it_object : it_scene.value().items()) {
							Object* tmp = nullptr;
							int objectType = it_object.value()["00 objectType"];
							if (objectType == 1) {
								vec3 posi;
								posi.x = it_object.value()["00 position"][0];
								posi.x = it_object.value()["00 position"][1];
								posi.x = it_object.value()["00 position"][2];
								double size;
								size = it_object.value()["01 size"];
								Material mat;
								mat.color.x = it_object.value()["02 material"]["color"][0];
								mat.color.y = it_object.value()["02 material"]["color"][1];
								mat.color.z = it_object.value()["02 material"]["color"][2];
								mat.emission.x = it_object.value()["02 material"]["emission"][0];
								mat.emission.y = it_object.value()["02 material"]["emission"][1];
								mat.emission.z = it_object.value()["02 material"]["emission"][2];
								tmp = new png::SphereObject(posi, size, mat);
							}
							data.object.push_back(tmp);
						}
					}
				}
			}
		}

	}

	Object::Object(const Material mat) {
		this->m_material = mat;
	}

	SphereObject::SphereObject(const vec3 position, const float size, const Material mat)
	: Object(mat)
	, m_position(position)
	, m_size(size)
	{}
	bool SphereObject::Intersect(const Ray& ray, double& out_dis, vec3 out_normal) const {
		const vec3 p_o = m_position - ray.org;
		const double b = Dot(p_o, ray.dir);
		const double D4 = b * b - Dot(p_o, p_o) + m_size * m_size;

		if (D4 < 0.0)
			return false;

		const double sqrt_D4 = sqrt(D4);
		const double t1 = b - sqrt_D4, t2 = b + sqrt_D4;

		const float minValue = 1e-5;
		if (t1 < minValue && t2 < minValue)
			return false;

		if (t1 > 0.001) {
			out_dis = t1;
			auto hitPoint = ray.org + ray.dir * out_dis;
			out_normal = Normalize(hitPoint - m_position);
			return true;
		}
		else {
			out_dis = t2;
			auto hitPoint = ray.org + ray.dir * out_dis;
			out_normal = Normalize(hitPoint - m_position);
			return true;
		}

		return false;
	}
	vec3 SphereObject::ComputeSurfacePoint(const std::function<double()> randGene) const
	{
		const auto theta = 2.0 * std::numbers::pi;
		const auto phi = 0.5 * std::numbers::pi;
		const vec3 localPoint = { std::sin(theta) * std::sin(phi), std::sin(theta) * std::cos(phi), std::cos(theta) };
		return localPoint * m_size + m_position;
	}

	PlaneObject::PlaneObject(const vec3 position, const vec3 up, const vec3 target, const double width, const Material mat)
	:Object(mat)
	,m_position(position)
	,m_up(up)
	{
		m_normal = Normalize(target - m_position);
		m_right = Normalize(Cross(m_normal, m_up)) * 0.5 * width;
	}
	bool PlaneObject::Intersect(const Ray& ray, double& out_dis, vec3 out_normal) const {
		if (std::abs(Dot(m_normal, ray.dir)) < FLT_EPSILON)
			return false;
	}
	vec3 PlaneObject::ComputeSurfacePoint(const std::function<double()> randGene) const {
		return vec3();
	}

}