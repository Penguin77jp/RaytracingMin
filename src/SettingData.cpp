#include "SettingData.h"
#include <fstream>
#include <iostream>
#include <format>

namespace png {
	LoadData::LoadData(std::string jsonName) {
		tinyxml2::XMLDocument doc;
		doc.LoadFile(jsonName.c_str());

		// innitialize strcut
		{
			data.renderType = -1;
			data.width = data.height = -1;
			data.samples = -1;
			data.spectrumSamples = -1;
			data.superSamples = -1;
			data.cameraOrigin = vec3();
			data.cameraTarget = vec3();
			data.fov = -1;
			data.object = std::vector<SceneObject*>();
		}

		from_xml(doc, data);
		Print();
	}

	void LoadData::to_xml(const SettingData& in_data, tinyxml2::XMLDocument& out_doc) {

	}
	void LoadData::from_xml(tinyxml2::XMLDocument& in_doc, SettingData& out_data) {
		in_doc.Print();
		std::function<std::optional<vec3>(tinyxml2::XMLElement*, std::string)> GetVec3FromXMLElement = [](tinyxml2::XMLElement* element, const std::string& str) {
			const char* val_str = nullptr;
			element->QueryStringAttribute(str.c_str(), &val_str);
			std::cout << "val_str" << val_str << std::endl;
			if (val_str == nullptr) {
				return std::optional<vec3>(std::nullopt);
			}
			else {
				vec3 ret;
				sscanf_s(val_str, "%lf,%lf,%lf", &ret.x, &ret.y, &ret.z);
				return std::optional<vec3>(ret);
			}
			};

		// render
		{
			const auto& render = in_doc.FirstChildElement("render");
			render->FindAttribute("type")->QueryIntValue(&out_data.renderType);
			render->FindAttribute("width")->QueryIntValue(&out_data.width);
			render->FindAttribute("height")->QueryIntValue(&out_data.height);
			{
				double tmp;
				render->FindAttribute("samples")->QueryDoubleValue(&tmp);
				out_data.samples = tmp;
			}
			render->FindAttribute("superSamples")->QueryIntValue(&out_data.superSamples);
		}
		// camera
		{
			const auto& camera = in_doc.FirstChildElement("camera");
			out_data.cameraOrigin = *GetVec3FromXMLElement(camera, "origin");
			out_data.cameraTarget = *GetVec3FromXMLElement(camera, "target");
			camera->FindAttribute("fov")->QueryDoubleValue(&out_data.fov);
		}
		// objects
		{
			tinyxml2::XMLElement* objs = in_doc.FirstChildElement("objects");
			tinyxml2::XMLElement* obj = objs->FirstChildElement("object");
			while (obj) {
				Material* mat = nullptr;
				{
					std::optional<vec3> color = GetVec3FromXMLElement(obj, "color");
					std::optional<vec3> emission_color = *GetVec3FromXMLElement(obj, "emission_color");
					if (!color.has_value() || !emission_color.has_value())
						break;
					mat = new DiffuseMaterial(new TextureSolid(*color), new TextureSolid(*emission_color), PDF_TYPE::CosImportSample);
				}
				const char* object_type = "";
				obj->QueryStringAttribute("type", &object_type);
				if (std::strcmp(object_type, "sphere") == 0) {
					std::optional<vec3> position = GetVec3FromXMLElement(obj, "position");
					double size = obj->FindAttribute("size")->DoubleValue();
					if (!position.has_value())
						break;
					SceneObject* tmp_object = new SphereObject(*position, size, mat);
					data.object.emplace_back(tmp_object);
				}
				else if (std::strcmp(object_type, "box")) {
					vec3 position = *GetVec3FromXMLElement(obj, "position");
					vec3 size = *GetVec3FromXMLElement(obj, "size");
					SceneObject* tmp_object = new BoxObject(position, size, mat);
				}
				else {
					std::cout << "invalid object type" << std::endl;
					break;
				}
				obj = obj->NextSiblingElement("object");
			}
		}

		/*
		for (auto& it : json.items()) {
			if (it.key() == "00 renderType") data.renderType = it.value();
			else if (it.key() == "00 width") data.width = it.value();
			else if (it.key() == "00 height") data.height = it.value();
			else if (it.key() == "00 samples") data.samples = it.value();
			else if (it.key() == "00 superSamples") data.superSamples = it.value();
			else if (it.key() == "00 spectrumSamples") data.spectrumSamples = it.value();
			else if (it.key() == "01 camera") {
				for (auto& it_cam : it.value().items()) {
					if (it_cam.key() == "origin") {
						data.cameraOrigin.x = it_cam.value()[0];
						data.cameraOrigin.y = it_cam.value()[1];
						data.cameraOrigin.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "target") {
						data.cameraTarget.x = it_cam.value()[0];
						data.cameraTarget.y = it_cam.value()[1];
						data.cameraTarget.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "fov") {
						data.fov = it_cam.value();
					}
				}
			}
			else if (it.key() == "02 scene") {
				for (auto& it_scene : it.value().items()) {
					if (it_scene.key() == "00 object") {
						for (auto& it_object : it_scene.value().items()) {
							//material
							int materialType = it_object.value()["03 material"]["type"];
							Material* mat = nullptr;
							if (materialType == 0) {
								vec3 color;
								color.x = it_object.value()["03 material"]["color"][0];
								color.y = it_object.value()["03 material"]["color"][1];
								color.z = it_object.value()["03 material"]["color"][2];
								TextureSolid* textureColor = new TextureSolid(color);
								vec3 emission;
								emission.x = it_object.value()["03 material"]["emission"][0];
								emission.y = it_object.value()["03 material"]["emission"][1];
								emission.z = it_object.value()["03 material"]["emission"][2];
								TextureSolid* textureEmission = new TextureSolid(emission);
								mat = new DiffuseMaterial(textureColor, textureEmission, PDF_TYPE::CosImportSample);
							}
							else if (materialType == 1) {
								vec3 color;
								color.x = it_object.value()["03 material"]["color"][0];
								color.y = it_object.value()["03 material"]["color"][1];
								color.z = it_object.value()["03 material"]["color"][2];
								TextureSolid* textureColor = new TextureSolid(color);
								vec3 emission;
								emission.x = it_object.value()["03 material"]["emission"][0];
								emission.y = it_object.value()["03 material"]["emission"][1];
								emission.z = it_object.value()["03 material"]["emission"][2];
								TextureSolid* textureEmission = new TextureSolid(emission);
								mat = new RefractionMaterial(TransparentMaterialType::BK7, textureColor, textureEmission);
							}

							SceneObject* tmp;
							int objectType = it_object.value()["00 objectType"];
							if (objectType == 0) {
								// sphere
								vec3 posi;
								posi.x = it_object.value()["01 position"][0];
								posi.y = it_object.value()["01 position"][1];
								posi.z = it_object.value()["01 position"][2];
								float size = it_object.value()["02 size"];
								tmp = new SphereObject(posi, size, mat);
							}
							else if (objectType == 1) {
								// mesh
								std::vector<vec3> posi;
								vec3 tmp_posi;
								int index = 0;
								for (auto& it_posi : it_object.value()["01 position"].items()) {
									if (index == 0) {
										tmp_posi.x = it_posi.value();
									}
									else if (index == 1) {
										tmp_posi.y = it_posi.value();
									}
									else if (index == 2) {
										tmp_posi.z = it_posi.value();
									}
									index++;
									if (index == 3) {
										posi.push_back(tmp_posi);
										index = 0;
									}
								}

								tmp = new MeshObject(posi, mat);
							}

							data.object.push_back(tmp);
						}
					}
				}
			}
		}
		*/

	}
	void LoadData::Print() {
		std::cout << "============== LoadData::Print() ==============" << std::endl;
		// render
		{
			std::cout << "render";
			std::cout << " type:" << data.renderType;
			std::cout << " width:" << data.width;
			std::cout << " height:" << data.height;
			std::cout << " samples:" << data.samples;
			std::cout << " superSamples:" << data.superSamples;
			std::cout << std::endl;
		}
		// camera
		{
			std::cout << "camera";
			std::cout << " origin:" << std::string(data.cameraOrigin);
			std::cout << " target:" << std::string(data.cameraTarget);
			std::cout << std::endl;
		}
		// objects
		{
			std::cout << "objects ";
			std::cout << "num:" << data.object.size() << std::endl;
			for (const auto& it : data.object) {
				std::cout << "object";
				SphereObject* sphere = dynamic_cast<SphereObject*>(it);
				if (sphere != nullptr) {
					std::cout << " sphere";
					std::cout << " posi:" << std::string(sphere->position());
					std::cout << " size:" << sphere->size();
					std::cout << std::endl;
				}
				else {
					std::cout << "invalid object" << std::endl;
				}

				std::cout << "material";
				DiffuseMaterial* diffuse = dynamic_cast<DiffuseMaterial*>(it->m_material);
				if (diffuse != nullptr) {
					std::cout << " diffuse";
					std::cout << " color:" << std::string(diffuse->m_color->GetColor(0, 0));
					std::cout << " emission:" << std::string(diffuse->m_emission->GetColor(0, 0));
					std::cout << std::endl;
				}
				else {
					std::cout << "invalid material" << std::endl;
				}
			}
		}

		std::cout << "============== END LoadData::Print() ==============" << std::endl;
	}
}