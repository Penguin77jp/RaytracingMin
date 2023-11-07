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
			data.samples = data.spectrumSamples = -1;
			data.superSamples = -1;
			data.cameraOrigin = vec3();
			data.cameraTarget = vec3();
			data.fov = -1;
		}

		from_xml(doc, data);
		Print();
	}

	void LoadData::to_xml(const SettingData& in_data, tinyxml2::XMLDocument& out_doc) {

	}
	void LoadData::from_xml(const tinyxml2::XMLDocument& in_doc, SettingData& out_data) {
		in_doc.Print();

		// render
		{
			const auto& render = in_doc.FirstChildElement("render");
			render->FindAttribute("type")->QueryIntValue(&out_data.renderType);
			render->FindAttribute("width")->QueryIntValue(&out_data.width);
			render->FindAttribute("height")->QueryIntValue(&out_data.height);
			render->FindAttribute("samples")->QueryIntValue(&out_data.samples);
			render->FindAttribute("superSamples")->QueryIntValue(&out_data.superSamples);
		}
		// camera
		{
			const auto& camera = in_doc.FirstChildElement("camera");
			const char* origin_xyz;
			camera->QueryStringAttribute("origin", &origin_xyz);
			sscanf_s(origin_xyz, "%lf,%lf,%lf", &out_data.cameraOrigin.x, &out_data.cameraOrigin.y, &out_data.cameraOrigin.z);
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
		std::cout << "============== END LoadData::Print() ==============" << std::endl;
	}
}